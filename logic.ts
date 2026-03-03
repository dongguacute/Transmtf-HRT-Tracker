// --- Types & Enums ---

export enum Route {
    injection = "injection",
    patchApply = "patchApply",
    patchRemove = "patchRemove",
    gel = "gel",
    oral = "oral",
    sublingual = "sublingual"
}

export enum Ester {
    E2 = "E2",
    EB = "EB",
    EV = "EV",
    EC = "EC",
    EN = "EN",
    CPA = "CPA"
}

export enum ExtraKey {
    concentrationMGmL = "concentrationMGmL",
    areaCM2 = "areaCM2",
    releaseRateUGPerDay = "releaseRateUGPerDay",
    sublingualTheta = "sublingualTheta",
    sublingualTier = "sublingualTier",
    gelSite = "gelSite"
}

enum GelSite {
    arm = "arm",
    thigh = "thigh",
    scrotal = "scrotal"
}

const GEL_SITE_ORDER = ["arm", "thigh", "scrotal"] as const;

const GelSiteParams = {
    [GelSite.arm]: 0.05,
    [GelSite.thigh]: 0.05,
    [GelSite.scrotal]: 0.40
};

export interface DoseEvent {
    id: string;
    route: Route;
    timeH: number; // Hours since 1970
    doseMG: number; // Dose in mg (of the ester/compound), NOT E2-equivalent
    ester: Ester;
    extras: Partial<Record<ExtraKey, number>>;
}

export interface SimulationResult {
    timeH: number[];
    concPGmL: number[];
    concPGmL_E2: number[];
    concPGmL_CPA: number[];
    auc: number;
}

// --- Lab Results & Calibration ---

export interface LabResult {
    id: string;
    timeH: number;
    concValue: number; // Value in the user's unit
    unit: 'pg/ml' | 'pmol/l';
}

export function convertToPgMl(val: number, unit: 'pg/ml' | 'pmol/l'): number {
    if (unit === 'pg/ml') return val;
    return val / 3.671; // pmol/L to pg/mL conversion
}

/**
 * Build a time-varying calibration scale based on lab results.
 * Returns a ratio function r(t) such that E2_conc(t) * r(t) is calibrated.
 * Strategy: compute ratio=obs/pred at each lab time, then linearly interpolate ratios over time.
 * NOTE: Lab results measure E2, not CPA, so calibration is only for E2.
 */
export function createCalibrationInterpolator(sim: SimulationResult | null, results: LabResult[]) {
    if (!sim || !results.length) return (_timeH: number) => 1;

    const getNearestConc_E2 = (timeH: number): number | null => {
        if (!sim.timeH.length) return null;
        let low = 0;
        let high = sim.timeH.length - 1;
        while (high - low > 1) {
            const mid = Math.floor((low + high) / 2);
            if (sim.timeH[mid] === timeH) return sim.concPGmL_E2[mid];
            if (sim.timeH[mid] < timeH) low = mid;
            else high = mid;
        }
        const idx = Math.abs(sim.timeH[high] - timeH) < Math.abs(sim.timeH[low] - timeH) ? high : low;
        return sim.concPGmL_E2[idx];
    };

    const points = results
        .map(r => {
            const obs = convertToPgMl(r.concValue, r.unit);
            let pred = interpolateConcentration_E2(sim, r.timeH);
            if (pred === null || Number.isNaN(pred)) {
                pred = getNearestConc_E2(r.timeH);
            }
            if (pred === null || pred <= 0.01 || obs <= 0) return null;
            const ratio = Math.max(0.1, Math.min(10, obs / pred));
            return { timeH: r.timeH, ratio };
        })
        .filter((p): p is { timeH: number; ratio: number } => !!p)
        .sort((a, b) => a.timeH - b.timeH);

    if (!points.length) return (_timeH: number) => 1;
    if (points.length === 1) {
        const r0 = points[0].ratio;
        return (_timeH: number) => r0;
    }

    return (timeH: number) => {
        if (timeH <= points[0].timeH) return points[0].ratio;
        if (timeH >= points[points.length - 1].timeH) return points[points.length - 1].ratio;
        // binary search
        let low = 0;
        let high = points.length - 1;
        while (high - low > 1) {
            const mid = Math.floor((low + high) / 2);
            if (points[mid].timeH === timeH) return points[mid].ratio;
            if (points[mid].timeH < timeH) low = mid;
            else high = mid;
        }
        const p1 = points[low];
        const p2 = points[high];
        const t = (timeH - p1.timeH) / (p2.timeH - p1.timeH);
        const r = p1.ratio + (p2.ratio - p1.ratio) * t;
        return Math.max(0.1, Math.min(10, r));
    };
}

/** Which Bayesian model to use for E2 calibration and CI bands. */
export type CalibrationModel = 'ekf' | 'ou-kalman';

// --- Bayesian OU-Kalman Calibration ---

/**
 * Parameters for the Ornstein-Uhlenbeck Kalman calibration model.
 *
 * The latent calibration log-factor x(t) follows an OU process:
 *   dx(t) = Θ(μ − x(t)) dt + σ dW_t
 *
 * The true E2 concentration is modelled as:
 *   c(t) = c₀(t) · exp(x(t))
 *
 * Lab observations in log-space:
 *   z_i = log(y_i) − log(c₀(t_i)) = x(t_i) + ε_i,   ε_i ~ N(0, τ²)
 */
export interface OUCalibParams {
    /** log-SD of measurement noise: sqrt(log(1 + CV²)). ~20% CV → 0.198. */
    tau: number;
    /** Mean-reversion rate (h⁻¹). Half-life = ln(2)/Theta. 7 days → ~0.00413 h⁻¹. */
    Theta: number;
    /** Process noise amplitude (log-scale h⁻⁰·⁵). Controls day-to-day variation. */
    sigma: number;
    /** Long-run log-mean (0 = unbiased, no systematic offset prior). */
    mu: number;
}

/** Default OU calibration parameters anchored to estradiol assay literature. */
export const OU_DEFAULT_PARAMS: OUCalibParams = {
    tau:   0.198,                     // ~20% CV assay noise (log-space SD)
    Theta: Math.LN2 / (7 * 24),       // 7-day mean-reversion half-life (~0.00413 h⁻¹)
    sigma: 0.02,                      // ~9% daily variation (log-scale)
    mu:    0.0,                       // unbiased prior
};

/**
 * Ornstein-Uhlenbeck Kalman filter + RTS smoother for Bayesian E2 calibration.
 *
 * Exact OU discretisation over interval dt:
 *   φ = exp(−Θ·dt)
 *   m⁻ = μ + φ·(m − μ)
 *   P⁻ = φ²·P + q,   q = P∞·(1 − φ²),   P∞ = σ²/(2Θ)
 *
 * Kalman update at observation z:
 *   K = P⁻ / (P⁻ + τ²)
 *   m = m⁻ + K·(z − m⁻)
 *   P = (1 − K)·P⁻
 *
 * RTS smoother (backward pass):
 *   G_i = P_i · φ_i / P⁻_{i+1}
 *   m_s[i] = m[i] + G_i · (m_s[i+1] − m⁻[i+1])
 *   P_s[i] = P[i] + G_i² · (P_s[i+1] − P⁻[i+1])
 *
 * @returns posterior mean m[] and variance P[] aligned with sim.timeH.
 */
export function buildOUKalmanCalibration(
    sim: SimulationResult,
    labResults: LabResult[],
    params: OUCalibParams = OU_DEFAULT_PARAMS
): { m: number[]; P: number[] } {
    const n = sim.timeH.length;
    if (n === 0) return { m: [], P: [] };

    const { tau, Theta, sigma, mu } = params;
    const tau2  = tau * tau;
    const P_inf = (sigma * sigma) / (2 * Theta);   // stationary variance of OU
    const EPS   = 0.1;                              // concentration floor (pg/mL)

    // Build log-space observations: z_i = log(y_i) − log(c₀(t_i))
    const tMin = sim.timeH[0];
    const tMax = sim.timeH[n - 1];
    const labs: { timeH: number; z: number }[] = [];
    for (const lab of labResults) {
        // Skip labs outside the simulation window — interpolation would clamp to
        // the edge value, producing a biased z that propagates via RTS smoothing.
        if (lab.timeH < tMin || lab.timeH > tMax) continue;
        const obs = convertToPgMl(lab.concValue, lab.unit);
        if (obs <= 0) continue;
        const c0 = interpolateConcentration_E2(sim, lab.timeH);
        if (c0 === null || c0 < EPS) continue;
        const z = Math.log(obs) - Math.log(c0);
        // Reject extreme outliers (> ±3.5 log-units ≈ ×33 ratio)
        if (!Number.isFinite(z) || Math.abs(z) > 3.5) continue;
        labs.push({ timeH: lab.timeH, z });
    }
    labs.sort((a, b) => a.timeH - b.timeH);

    // Build combined time grid: sim.timeH ∪ {lab.timeH} — sorted unique
    const gridSet = new Set<number>(sim.timeH);
    for (const lab of labs) gridSet.add(lab.timeH);
    const grid = Array.from(gridSet).sort((a, b) => a - b);
    const ng = grid.length;

    // O(1) lookup from time → grid index
    const gridIndex = new Map<number, number>();
    for (let i = 0; i < ng; i++) gridIndex.set(grid[i], i);

    // ── Forward Kalman filter ──────────────────────────────────────────
    const mFwd  = new Float64Array(ng).fill(mu);
    const PFwd  = new Float64Array(ng).fill(P_inf);
    const mPred = new Float64Array(ng).fill(mu);
    const PPred = new Float64Array(ng).fill(P_inf);

    let m = mu;
    let P = P_inf;
    let labPtr = 0;

    for (let i = 0; i < ng; i++) {
        // Predict (OU transition from previous time)
        if (i > 0) {
            const dt = grid[i] - grid[i - 1];
            if (dt > 0) {
                const phi = Math.exp(-Theta * dt);
                const q   = P_inf * (1 - phi * phi);
                m = mu + phi * (m - mu);
                P = phi * phi * P + q;
            }
        }
        mPred[i] = m;
        PPred[i] = P;

        // Update with any lab(s) at exactly this grid time
        while (labPtr < labs.length && labs[labPtr].timeH === grid[i]) {
            const S = P + tau2;
            const K = P / S;
            m = m + K * (labs[labPtr].z - m);
            P = (1 - K) * P;
            labPtr++;
        }
        mFwd[i] = m;
        PFwd[i] = Math.max(P, 1e-12);
    }

    // ── RTS smoother (backward pass) ──────────────────────────────────
    const mS = Float64Array.from(mFwd);
    const PS = Float64Array.from(PFwd);

    for (let i = ng - 2; i >= 0; i--) {
        const dt = grid[i + 1] - grid[i];
        if (dt <= 0) continue;
        const phi = Math.exp(-Theta * dt);
        // G = P_i · φ / P⁻_{i+1}  (Kalman smoother gain)
        const G = PPred[i + 1] > 1e-12 ? PFwd[i] * phi / PPred[i + 1] : 0;
        mS[i] = mFwd[i] + G * (mS[i + 1] - mPred[i + 1]);
        PS[i] = Math.max(PFwd[i] + G * G * (PS[i + 1] - PPred[i + 1]), 1e-9);
    }

    // ── Extract at sim.timeH positions (all sim times are in gridSet) ─
    const outM = new Array<number>(n);
    const outP = new Array<number>(n);
    for (let j = 0; j < n; j++) {
        const idx = gridIndex.get(sim.timeH[j]);
        if (idx !== undefined) {
            outM[j] = mS[idx];
            outP[j] = PS[idx];
        } else {
            // Defensive fallback (should never occur)
            outM[j] = mu;
            outP[j] = P_inf;
        }
    }

    return { m: outM, P: outP };
}

// --- Constants & Parameters (PKparameter.swift & PKcore.swift) ---

const CorePK = {
    vdPerKG: 2.0, // L/kg for E2
    /** @deprecated Use CPA_2COMP_PK.V1_per_kg (central Vc) instead of this apparent Vd */
    vdPerKG_CPA: 14.0, // L/kg apparent Vd (= Vss ≈ 986 L/70 kg) — NOT used for concentration
    kClear: 0.41,
    kClearInjection: 0.041,
    depotK1Corr: 1.0
};

/**
 * CPA 2-compartment oral PK constants derived from SmPC (Diane-35 / Androcur).
 *
 * Macro constants α, β and micro constants are derived from:
 *   CL  = 3.6 mL/min/kg  → 0.216 L/h/kg  (≈15.12 L/h for 70 kg)
 *   Vss = 986 L / 70 kg  ≈ 14.09 L/kg
 *   α   = 0.866 h⁻¹      (t1/2 ≈ 0.8 h, fast distribution phase)
 *   β   = 0.01031 h⁻¹    (t1/2 ≈ 67 h ≈ 2.8 days, terminal elimination)
 *   k21 = α·β / (α+β − CL/Vc)    where Vc = CL/k10
 *   V1  = CL / k10 = 0.734 L/kg  (central compartment, ≈51 L for 70 kg)
 *
 * References: EMA/Bayer SmPC for cyproterone-containing products.
 */
const CPA_2COMP_PK = {
    F: 0.88,          // oral bioavailability (~88% per SmPC; old code used 0.70)
    ka: 0.5,          // oral absorption rate h⁻¹  →  Tmax ≈ 1.5 h
    alpha: 0.8660,    // fast hybrid rate constant h⁻¹  (t1/2 ≈ 0.8 h)
    beta: 0.01031,    // slow hybrid rate constant h⁻¹  (t1/2 ≈ 67 h ≈ 2.8 days)
    k21: 0.03034,     // peripheral→central transfer h⁻¹
    V1_per_kg: 0.734, // central compartment Vc per kg  (L/kg);  Vc ≈ 51 L for 70 kg
    popLogVar: 0.09,  // population PK log-variance ≈ (0.30)² — ~30% CV for Vd/CL uncertainty
};

const EsterInfo = {
    [Ester.E2]: { name: "Estradiol", mw: 272.38 },
    [Ester.EB]: { name: "Estradiol Benzoate", mw: 376.50 },
    [Ester.EV]: { name: "Estradiol Valerate", mw: 356.50 },
    [Ester.EC]: { name: "Estradiol Cypionate", mw: 396.58 },
    [Ester.EN]: { name: "Estradiol Enanthate", mw: 384.56 },
    [Ester.CPA]: { name: "Cyproterone Acetate", mw: 416.94 }
};

export function getToE2Factor(ester: Ester): number {
    if (ester === Ester.E2) return 1.0;
    return EsterInfo[Ester.E2].mw / EsterInfo[ester].mw;
}

const TwoPartDepotPK = {
    Frac_fast: { [Ester.EB]: 0.90, [Ester.EV]: 0.40, [Ester.EC]: 0.229164549, [Ester.EN]: 0.05, [Ester.E2]: 1.0 },
    k1_fast: { [Ester.EB]: 0.144, [Ester.EV]: 0.0216, [Ester.EC]: 0.005035046, [Ester.EN]: 0.0010, [Ester.E2]: 0.5 }, // Added non-zero k1 for E2
    k1_slow: { [Ester.EB]: 0.114, [Ester.EV]: 0.0138, [Ester.EC]: 0.004510574, [Ester.EN]: 0.0050, [Ester.E2]: 0 }
};

const InjectionPK = {
    formationFraction: { [Ester.EB]: 0.1092, [Ester.EV]: 0.0623, [Ester.EC]: 0.1173, [Ester.EN]: 0.12, [Ester.E2]: 1.0 }
};

const EsterPK = {
    k2: { [Ester.EB]: 0.090, [Ester.EV]: 0.070, [Ester.EC]: 0.045, [Ester.EN]: 0.015, [Ester.E2]: 0 }
};

const OralPK = {
    kAbsE2: 0.32,
    kAbsEV: 0.05,
    bioavailability: 0.03,
    kAbsSL: 1.8
};

// Define deterministic order for mapping integer tiers (0-3)   to keys
export const SL_TIER_ORDER = ["quick", "casual", "standard", "strict"] as const;

export const SublingualTierParams = {
    quick: { theta: 0.01, hold: 2 },
    casual: { theta: 0.04, hold: 5 },
    standard: { theta: 0.11, hold: 10 },
    strict: { theta: 0.18, hold: 15 }
};

export function getBioavailabilityMultiplier(
    route: Route,
    ester: Ester,
    extras: Partial<Record<ExtraKey, number>> = {}
): number {
    const mwFactor = getToE2Factor(ester);
    
    switch (route) {
        case Route.injection: {
            const formation = InjectionPK.formationFraction[ester] ?? 0.08;
            return formation * mwFactor;
        }
        case Route.oral:
            return OralPK.bioavailability * mwFactor;
        case Route.sublingual: {
            let theta = 0.11;
            if (extras[ExtraKey.sublingualTheta] !== undefined) {
                const customTheta = extras[ExtraKey.sublingualTheta];
                if (typeof customTheta === 'number' && Number.isFinite(customTheta)) {
                    theta = Math.min(1, Math.max(0, customTheta));
                }
            } else if (extras[ExtraKey.sublingualTier] !== undefined) {
                const tierIdx = Math.min(SL_TIER_ORDER.length - 1, Math.max(0, Math.round(extras[ExtraKey.sublingualTier]!)));
                const tierKey = SL_TIER_ORDER[tierIdx] || 'standard';
                theta = SublingualTierParams[tierKey]?.theta ?? 0.11;
            }
            return (theta + (1 - theta) * OralPK.bioavailability) * mwFactor;
        }
        case Route.gel: {
            const siteIdx = Math.min(GEL_SITE_ORDER.length - 1, Math.max(0, Math.round(extras[ExtraKey.gelSite] ?? 0)));
            // @ts-ignore
            const siteKey = GEL_SITE_ORDER[siteIdx] || GelSite.arm;
            const bio = GelSiteParams[siteKey] ?? 0.05;
            return bio * mwFactor;
        }
        case Route.patchApply:
            return 1.0 * mwFactor;
        case Route.patchRemove:
        default:
            return 0;
    }
}

// --- Math Models ---

interface PKParams {
    Frac_fast: number;
    k1_fast: number;
    k1_slow: number;
    k2: number;
    k3: number;
    F: number;
    rateMGh: number;
    F_fast: number;
    F_slow: number;
}

function resolveParams(event: DoseEvent): PKParams {
    const defaultK3 = event.route === Route.injection ? CorePK.kClearInjection : CorePK.kClear;
    const toE2 = getToE2Factor(event.ester);
    const extras = event.extras ?? {};

    switch (event.route) {
        case Route.injection: {
            const Frac_fast = TwoPartDepotPK.Frac_fast[event.ester] ?? 0.5;
            const k1_fast = (TwoPartDepotPK.k1_fast[event.ester] ?? 0.1) * CorePK.depotK1Corr;
            const k1_slow = (TwoPartDepotPK.k1_slow[event.ester] ?? 0.01) * CorePK.depotK1Corr;
            const k2 = EsterPK.k2[event.ester] ?? 0;
            const F = getBioavailabilityMultiplier(Route.injection, event.ester, extras);
            return { Frac_fast, k1_fast, k1_slow, k2, k3: defaultK3, F, rateMGh: 0, F_fast: F, F_slow: F };
        }

        case Route.sublingual: {
            let theta = 0.11;
            if (extras[ExtraKey.sublingualTheta] !== undefined) {
                const customTheta = extras[ExtraKey.sublingualTheta];
                if (typeof customTheta === 'number' && Number.isFinite(customTheta)) {
                    theta = Math.min(1, Math.max(0, customTheta));
                }
            } else if (extras[ExtraKey.sublingualTier] !== undefined) {
                const tierRaw = extras[ExtraKey.sublingualTier];
                if (typeof tierRaw === 'number' && Number.isFinite(tierRaw)) {
                    const tierIdx = Math.min(SL_TIER_ORDER.length - 1, Math.max(0, Math.round(tierRaw)));
                    const tierKey = SL_TIER_ORDER[tierIdx] || 'standard';
                    theta = SublingualTierParams[tierKey]?.theta ?? theta;
                }
            }
            const k1_fast = OralPK.kAbsSL;
            const k1_slow = event.ester === Ester.EV ? OralPK.kAbsEV : OralPK.kAbsE2;
            const k2 = EsterPK.k2[event.ester] ?? 0;
            const F_fast = toE2;
            const F_slow = OralPK.bioavailability * toE2;
            const F = theta * F_fast + (1 - theta) * F_slow;
            return { Frac_fast: theta, k1_fast, k1_slow, k2, k3: defaultK3, F, rateMGh: 0, F_fast, F_slow };
        }

        case Route.gel: {
            const F = getBioavailabilityMultiplier(Route.gel, event.ester, extras);
            const k1 = 0.022;
            return { Frac_fast: 1.0, k1_fast: k1, k1_slow: 0, k2: 0, k3: defaultK3, F, rateMGh: 0, F_fast: F, F_slow: F };
        }

        case Route.patchApply: {
            const F = getBioavailabilityMultiplier(Route.patchApply, event.ester, extras);
            const releaseRateUGPerDay = extras[ExtraKey.releaseRateUGPerDay];
            const rateMGh = (typeof releaseRateUGPerDay === 'number' && Number.isFinite(releaseRateUGPerDay) && releaseRateUGPerDay > 0)
                ? (releaseRateUGPerDay / 24 / 1000) * F
                : 0;
            if (rateMGh > 0) {
                return { Frac_fast: 1.0, k1_fast: 0, k1_slow: 0, k2: 0, k3: defaultK3, F, rateMGh, F_fast: F, F_slow: F };
            }
            const k1 = 0.0075;
            return { Frac_fast: 1.0, k1_fast: k1, k1_slow: 0, k2: 0, k3: defaultK3, F, rateMGh: 0, F_fast: F, F_slow: F };
        }

        case Route.patchRemove:
            return { Frac_fast: 0, k1_fast: 0, k1_slow: 0, k2: 0, k3: defaultK3, F: 0, rateMGh: 0, F_fast: 0, F_slow: 0 };

        case Route.oral: {
            // === 针对 CPA 的特殊处理开始 ===
            if (event.ester === Ester.CPA) {
                return {
                    Frac_fast: 1.0,
                    k1_fast: 1.0,
                    k1_slow: 0,
                    k2: 0,
                    k3: 0.017,
                    F: 0.7,
                    rateMGh: 0,
                    F_fast: 0.7,
                    F_slow: 0.7
                };
            }
            // === 针对 CPA 的特殊处理结束 ===

            const k1Value = event.ester === Ester.EV ? OralPK.kAbsEV : OralPK.kAbsE2;
            const k2Value = event.ester === Ester.EV ? (EsterPK.k2[Ester.EV] || 0) : 0;
            const F = OralPK.bioavailability * toE2;
            return { Frac_fast: 1.0, k1_fast: k1Value, k1_slow: 0, k2: k2Value, k3: defaultK3, F, rateMGh: 0, F_fast: F, F_slow: F };
        }
    }

    return { Frac_fast: 0, k1_fast: 0, k1_slow: 0, k2: 0, k3: defaultK3, F: 0, rateMGh: 0, F_fast: 0, F_slow: 0 };
}

/**
 * Compute CPA amount in the central compartment (mg) at time tau after a single
 * oral dose using the 2-compartment model derived from SmPC PK parameters.
 *
 * Analytical solution for central compartment amount (X₁) with first-order oral absorption:
 *   X₁(t) = F·D·ka · [
 *     (k21−ka) / ((α−ka)(β−ka)) · exp(−ka·t)  +
 *     (k21−α)  / ((ka−α)(β−α)) · exp(−α·t)   +
 *     (k21−β)  / ((ka−β)(α−β)) · exp(−β·t)
 *   ]
 *
 * Dividing by V1 gives central-compartment concentration.
 */
function compute2CompCPACentralAmount(doseMG: number, tau: number): number {
    if (tau < 0 || doseMG <= 0) return 0;
    const { F, ka, alpha, beta, k21 } = CPA_2COMP_PK;
    // Guard against near-singularity (should not occur with fixed SmPC params, but defensive)
    const eps = 1e-8;
    if (Math.abs(alpha - ka) < eps || Math.abs(beta - ka) < eps || Math.abs(alpha - beta) < eps) {
        // Fallback: single-compartment with terminal rate (conservative approximation)
        if (Math.abs(ka - beta) < eps) return Math.max(0, doseMG * F * ka * tau * Math.exp(-beta * tau));
        return Math.max(0, doseMG * F * ka / (ka - beta) * (Math.exp(-beta * tau) - Math.exp(-ka * tau)));
    }
    const A = (k21 - ka)    / ((alpha - ka) * (beta - ka));
    const B = (k21 - alpha) / ((ka - alpha) * (beta - alpha));
    const C = (k21 - beta)  / ((ka - beta)  * (alpha - beta));
    const val = doseMG * F * ka * (
        A * Math.exp(-ka    * tau) +
        B * Math.exp(-alpha * tau) +
        C * Math.exp(-beta  * tau)
    );
    return Math.max(0, val);
}

// 3-Compartment Analytical Solution
function _analytic3C(tau: number, doseMG: number, F: number, k1: number, k2: number, k3: number): number {
    if (k1 <= 0 || doseMG <= 0) return 0;
    const k1_k2 = k1 - k2;
    const k1_k3 = k1 - k3;
    const k2_k3 = k2 - k3;

    if (Math.abs(k1_k2) < 1e-9 || Math.abs(k1_k3) < 1e-9 || Math.abs(k2_k3) < 1e-9) return 0; // Singularity protection

    const term1 = Math.exp(-k1 * tau) / (k1_k2 * k1_k3);
    const term2 = Math.exp(-k2 * tau) / (-k1_k2 * k2_k3);
    const term3 = Math.exp(-k3 * tau) / (k1_k3 * k2_k3);

    return doseMG * F * k1 * k2 * (term1 + term2 + term3);
}

function oneCompAmount(tau: number, doseMG: number, p: PKParams): number {
    const k1 = p.k1_fast;
    if (Math.abs(k1 - p.k3) < 1e-9) {
        return doseMG * p.F * k1 * tau * Math.exp(-p.k3 * tau);
    }
    return doseMG * p.F * k1 / (k1 - p.k3) * (Math.exp(-p.k3 * tau) - Math.exp(-k1 * tau));
}

// Model Solver
class PrecomputedEventModel {
    private model: (t: number) => number;

    constructor(event: DoseEvent, allEvents: DoseEvent[]) {
        const params = resolveParams(event);
        const startTime = event.timeH;
        const dose = event.doseMG;

        switch (event.route) {
            case Route.injection:
                this.model = (timeH: number) => {
                    const tau = timeH - startTime;
                    if (tau < 0) return 0;
                         const doseFast = dose * params.Frac_fast;
                         const doseSlow = dose * (1.0 - params.Frac_fast);

                         return _analytic3C(tau, doseFast, params.F, params.k1_fast, params.k2, params.k3) +
                             _analytic3C(tau, doseSlow, params.F, params.k1_slow, params.k2, params.k3);
                };
                break;
            case Route.gel:
            case Route.oral:
                this.model = (timeH: number) => {
                    const tau = timeH - startTime;
                    if (tau < 0) return 0;
                    // CPA oral uses dedicated 2-compartment model
                    if (event.ester === Ester.CPA) {
                        return compute2CompCPACentralAmount(dose, tau);
                    }
                    return oneCompAmount(tau, dose, params);
                };
                break;
            case Route.sublingual:
                this.model = (timeH: number) => {
                    const tau = timeH - startTime;
                    if (tau < 0) return 0;
                    if (params.k2 > 0) {
                        // EV Sublingual
                        const doseF = dose * params.Frac_fast;
                        const doseS = dose * (1.0 - params.Frac_fast);
                        return _analytic3C(tau, doseF, params.F_fast, params.k1_fast, params.k2, params.k3) +
                               _analytic3C(tau, doseS, params.F_slow, params.k1_slow, params.k2, params.k3);
                    } else {
                        // E2 Sublingual
                        const doseF = dose * params.Frac_fast;
                        const doseS = dose * (1.0 - params.Frac_fast);
                        
                        // Helper for dual branch 1st order
                        const branch = (d: number, F: number, ka: number, ke: number, t: number) => {
                             if (Math.abs(ka - ke) < 1e-9) return d * F * ka * t * Math.exp(-ke * t);
                             return d * F * ka / (ka - ke) * (Math.exp(-ke * t) - Math.exp(-ka * t));
                        };
                        return branch(doseF, params.F_fast, params.k1_fast, params.k3, tau) +
                               branch(doseS, params.F_slow, params.k1_slow, params.k3, tau);
                    }
                };
                break;
            case Route.patchApply:
                const remove = allEvents.find(e => e.route === Route.patchRemove && e.timeH > startTime);
                const wearH = (remove?.timeH ?? Number.MAX_VALUE) - startTime;
                
                this.model = (timeH: number) => {
                    const tau = timeH - startTime;
                    if (tau < 0) return 0;
                    
                    // Zero Order
                    if (params.rateMGh > 0) {
                        if (tau <= wearH) {
                            return params.rateMGh / params.k3 * (1 - Math.exp(-params.k3 * tau));
                        } else {
                            const amtRemoval = params.rateMGh / params.k3 * (1 - Math.exp(-params.k3 * wearH));
                            return amtRemoval * Math.exp(-params.k3 * (tau - wearH));
                        }
                    }
                    // First order legacy
                    const amtUnderPatch = oneCompAmount(tau, dose, params);
                    if (tau > wearH) {
                        const amtAtRemoval = oneCompAmount(wearH, dose, params);
                        return amtAtRemoval * Math.exp(-params.k3 * (tau - wearH));
                    }
                    return amtUnderPatch;
                };
                break;
            default:
                this.model = () => 0;
        }
    }

    amount(timeH: number): number {
        return this.model(timeH);
    }
}

// --- Simulation Engine ---

export function runSimulation(events: DoseEvent[], bodyWeightKG: number): SimulationResult | null {
    if (events.length === 0) return null;

    const sortedEvents = [...events].sort((a, b) => a.timeH - b.timeH);
    const precomputed = sortedEvents
        .filter(e => e.route !== Route.patchRemove)
        .map(e => ({ model: new PrecomputedEventModel(e, sortedEvents), ester: e.ester }));

    const startTime = sortedEvents[0].timeH - 24;
    const endTime = sortedEvents[sortedEvents.length - 1].timeH + (24 * 14);
    const steps = 1000;

    // Different Vd for E2 and CPA
    const plasmaVolumeML_E2 = CorePK.vdPerKG * bodyWeightKG * 1000; // E2: ~2.0 L/kg
    // CPA: use central compartment V1, not apparent Vss — 2-compartment model requires this
    const plasmaVolumeML_CPA = CPA_2COMP_PK.V1_per_kg * bodyWeightKG * 1000; // Vc ≈ 51 L for 70 kg

    const timeH: number[] = [];
    const concPGmL: number[] = [];
    const concPGmL_E2: number[] = [];
    const concPGmL_CPA: number[] = []; // Will store in ng/mL (not pg/mL)
    let auc = 0;

    const stepSize = (endTime - startTime) / (steps - 1);

    for (let i = 0; i < steps; i++) {
        const t = startTime + i * stepSize;
        let totalAmountMG_E2 = 0;
        let totalAmountMG_CPA = 0;

        for (const { model, ester } of precomputed) {
            const amount = model.amount(t);
            if (ester === Ester.CPA) {
                totalAmountMG_CPA += amount;
            } else {
                totalAmountMG_E2 += amount;
            }
        }

        // E2: pg/mL (using E2 Vd)
        const currentConc_E2 = (totalAmountMG_E2 * 1e9) / plasmaVolumeML_E2;

        // CPA: ng/mL (using CPA Vd, convert from mg to ng: 1e6 instead of 1e9)
        const currentConc_CPA = (totalAmountMG_CPA * 1e6) / plasmaVolumeML_CPA;

        // Total in pg/mL (convert CPA from ng/mL to pg/mL for compatibility)
        const currentConc = currentConc_E2 + (currentConc_CPA * 1000);

        timeH.push(t);
        concPGmL.push(currentConc);
        concPGmL_E2.push(currentConc_E2); // pg/mL
        concPGmL_CPA.push(currentConc_CPA); // ng/mL

        if (i > 0) {
            auc += 0.5 * (currentConc + concPGmL[i - 1]) * stepSize;
        }
    }

    return { timeH, concPGmL, concPGmL_E2, concPGmL_CPA, auc };
}

export function interpolateConcentration(sim: SimulationResult, hour: number): number | null {
    if (!sim.timeH.length) return null;
    if (hour <= sim.timeH[0]) return sim.concPGmL[0];
    if (hour >= sim.timeH[sim.timeH.length - 1]) return sim.concPGmL[sim.concPGmL.length - 1];

    // Binary search for efficiency
    let low = 0;
    let high = sim.timeH.length - 1;

    while (high - low > 1) {
        const mid = Math.floor((low + high) / 2);
        if (sim.timeH[mid] === hour) return sim.concPGmL[mid];
        if (sim.timeH[mid] < hour) low = mid;
        else high = mid;
    }

    const t0 = sim.timeH[low];
    const t1 = sim.timeH[high];
    const c0 = sim.concPGmL[low];
    const c1 = sim.concPGmL[high];

    if (t1 === t0) return c0;
    const ratio = (hour - t0) / (t1 - t0);
    return c0 + (c1 - c0) * ratio;
}

export function interpolateConcentration_E2(sim: SimulationResult, hour: number): number | null {
    if (!sim.timeH.length) return null;
    if (hour <= sim.timeH[0]) return sim.concPGmL_E2[0];
    if (hour >= sim.timeH[sim.timeH.length - 1]) return sim.concPGmL_E2[sim.concPGmL_E2.length - 1];

    // Binary search for efficiency
    let low = 0;
    let high = sim.timeH.length - 1;

    while (high - low > 1) {
        const mid = Math.floor((low + high) / 2);
        if (sim.timeH[mid] === hour) return sim.concPGmL_E2[mid];
        if (sim.timeH[mid] < hour) low = mid;
        else high = mid;
    }

    const t0 = sim.timeH[low];
    const t1 = sim.timeH[high];
    const c0 = sim.concPGmL_E2[low];
    const c1 = sim.concPGmL_E2[high];

    if (t1 === t0) return c0;
    const ratio = (hour - t0) / (t1 - t0);
    return c0 + (c1 - c0) * ratio;
}

export function interpolateConcentration_CPA(sim: SimulationResult, hour: number): number | null {
    if (!sim.timeH.length) return null;
    if (hour <= sim.timeH[0]) return sim.concPGmL_CPA[0];
    if (hour >= sim.timeH[sim.timeH.length - 1]) return sim.concPGmL_CPA[sim.concPGmL_CPA.length - 1];

    // Binary search for efficiency
    let low = 0;
    let high = sim.timeH.length - 1;

    while (high - low > 1) {
        const mid = Math.floor((low + high) / 2);
        if (sim.timeH[mid] === hour) return sim.concPGmL_CPA[mid];
        if (sim.timeH[mid] < hour) low = mid;
        else high = mid;
    }

    const t0 = sim.timeH[low];
    const t1 = sim.timeH[high];
    const c0 = sim.concPGmL_CPA[low];
    const c1 = sim.concPGmL_CPA[high];

    if (t1 === t0) return c0;
    const ratio = (hour - t0) / (t1 - t0);
    return c0 + (c1 - c0) * ratio;
}

// --- Personal PK Learning System (EKF / MAP-Bayesian) ---
//
// Implements a 2-parameter Extended Kalman Filter (EKF) that learns individual
// PK parameters from lab calibration points.
//
// State vector theta = [theta_s, theta_k]:
//   theta_s → amplitude scale  s = exp(theta_s)   (adjusts Vd / systemic exposure)
//   theta_k → clearance scale k = exp(theta_k)   (adjusts kClear & kClearInjection)
//
// Prediction in log-space: yhat(t) = log(s * C_pk(t, k)) = theta_s + log(C_pk(t,k))
// This makes the gradient ∂yhat/∂theta_s = 1 analytically, reducing finite-diff calls.

export interface ResidualAnchor {
    timeH: number;
    logRatio: number; // log(observed_pgml) - log(predicted_with_theta_pgml)
    w: number;        // confidence weight [0, 1]
    kind: 'lab';
}

export interface PersonalModelState {
    modelVersion: 'pk-ekf-v1';
    thetaMean: [number, number];                       // [theta_s, theta_k]
    thetaCov: [[number, number], [number, number]];    // 2×2 posterior covariance
    Q: [[number, number], [number, number]];           // process noise (slow drift)
    Rlog: number;                                      // measurement noise (log-space)
    anchors: ResidualAnchor[];
    observationCount: number;
    updatedAt: string;                                 // ISO-8601
}

export interface EKFDiagnostics {
    NIS: number;              // Normalized Innovation Squared (chi-sq test stat)
    isOutlier: boolean;       // NIS > chi2 at 95% (1 DOF) → flag this observation
    residualLog: number;      // y - yhat in log-space (positive = obs > pred)
    predictedPGmL: number;   // E2 prediction before update
    observedPGmL: number;    // E2 observation
    ci95Low: number;          // 95% CI lower bound (pg/mL) at observation time
    ci95High: number;         // 95% CI upper bound (pg/mL) at observation time
    convergenceScore: number; // 0 = no convergence, 1 = fully converged
    thetaS: number;           // exp(theta_s): current amplitude multiplier
    thetaK: number;           // exp(theta_k): current clearance scale
}

// Default prior parameters (see paper: sigma_s=0.5, sigma_k=0.3, sigma_y≈0.2)
const EKF_INITIAL_COV: [[number, number], [number, number]] = [[0.25, 0.0], [0.0, 0.09]];
const EKF_Q: [[number, number], [number, number]] = [[0.0004, 0.0], [0.0, 0.0001]];
const EKF_RLOG = 0.04;                    // (0.2 log-space SD)²
const EKF_EPS = 0.1;                      // concentration floor (pg/mL)
const EKF_EPS_CPA = 0.001;                // concentration floor for CPA (ng/mL)
const EKF_CHI2_95 = 3.841;               // chi-squared 95th percentile, 1 DOF
const EKF_DELTA_K = 0.01;                // finite-difference step for theta_k
const EKF_CI_MAX_E2 = 5000;               // chart safety cap (pg/mL)
const EKF_CI_MAX_CPA = 500;               // chart safety cap (ng/mL)
/** Irreducible predictive log-SD (assay CV ~5–15% + biological variability ~10–20% + model error).
 *  σ = 0.27 → 95% predictive band ≈ ×/÷ 1.70.  Range 0.25–0.32 is defensible. */
const EKF_SIGMA_RESIDUAL_LOG = 0.27;
/** Reference period for time-scaled Q drift (30 days). EKF_Q is calibrated for this interval. */
const EKF_Q_REF_PERIOD_H = 30 * 24;

/**
 * CPA → E2 clearance inhibition model parameters (Chou–Talalay / Hill).
 *
 * Scientific basis: CPA is metabolised by CYP3A4; competitive saturation of
 * CYP3A4 may reduce E2 first-pass clearance when CPA concentrations are high.
 * Parameters are conservative estimates — clinical magnitude is debated.
 * Feature is off by default; enabled via applyCPAInhibitionToE2.
 *
 * Note on saturation: at typical HRT doses (10–50 mg/day oral CPA), steady-state
 * Cavg ≈ 24–121 ng/mL (using 2-comp model). With IC50=50 ng/mL and Emax=0.20,
 * inhibition ranges from ~0.09 to ~0.17 — avoiding premature saturation.
 * IC50 was deliberately set above the typical HRT trough (~5–15 ng/mL) to keep
 * the model in the dynamic (non-saturated) range for clinically relevant doses.
 */
const CPA_E2_INHIBITION = {
    Emax: 0.20,   // maximum fraction of E2 clearance inhibited (20%)
    IC50: 50.0,   // CPA ng/mL for 50% of Emax — above typical HRT trough to stay in dynamic range
    n:    1.0,    // Hill coefficient (1 = Michaelis-Menten)
};

/** Initialise a new personal model state at the population prior. */
export function initPersonalModel(): PersonalModelState {
    return {
        modelVersion: 'pk-ekf-v1',
        thetaMean: [0, 0],
        thetaCov: EKF_INITIAL_COV,
        Q: EKF_Q,
        Rlog: EKF_RLOG,
        anchors: [],
        observationCount: 0,
        updatedAt: new Date().toISOString(),
    };
}

/**
 * Compute E2 amount contributed by one event at time tau after the event,
 * using a scaled clearance rate (kScale = exp(theta_k)).
 * Only E2-family compounds are evaluated; CPA events return 0.
 */
function computeEventAmountWithKScale(
    event: DoseEvent,
    allEvents: DoseEvent[],
    tau: number,
    kScale: number
): number {
    if (tau < 0) return 0;
    if (event.route === Route.patchRemove) return 0;
    if (event.ester === Ester.CPA) return 0;

    const params = resolveParams(event);
    const k3 = params.k3 * kScale; // scale the elimination rate

    switch (event.route) {
        case Route.injection: {
            const doseFast = event.doseMG * params.Frac_fast;
            const doseSlow = event.doseMG * (1.0 - params.Frac_fast);
            return _analytic3C(tau, doseFast, params.F, params.k1_fast, params.k2, k3) +
                   _analytic3C(tau, doseSlow, params.F, params.k1_slow, params.k2, k3);
        }
        case Route.gel:
        case Route.oral: {
            const paramsK = { ...params, k3 };
            return oneCompAmount(tau, event.doseMG, paramsK);
        }
        case Route.sublingual: {
            const doseF = event.doseMG * params.Frac_fast;
            const doseS = event.doseMG * (1.0 - params.Frac_fast);
            if (params.k2 > 0) {
                // EV sublingual (3-compartment each branch)
                return _analytic3C(tau, doseF, params.F_fast, params.k1_fast, params.k2, k3) +
                       _analytic3C(tau, doseS, params.F_slow, params.k1_slow, params.k2, k3);
            } else {
                // E2 sublingual (1-compartment dual branch)
                const branch = (d: number, F: number, ka: number, ke: number, t: number): number => {
                    if (Math.abs(ka - ke) < 1e-9) return d * F * ka * t * Math.exp(-ke * t);
                    return d * F * ka / (ka - ke) * (Math.exp(-ke * t) - Math.exp(-ka * t));
                };
                return branch(doseF, params.F_fast, params.k1_fast, k3, tau) +
                       branch(doseS, params.F_slow, params.k1_slow, k3, tau);
            }
        }
        case Route.patchApply: {
            const remove = allEvents.find(e => e.route === Route.patchRemove && e.timeH > event.timeH);
            const wearH = (remove?.timeH ?? Number.MAX_VALUE) - event.timeH;
            if (params.rateMGh > 0) {
                if (tau <= wearH) {
                    return params.rateMGh / k3 * (1 - Math.exp(-k3 * tau));
                } else {
                    const amtAtRemoval = params.rateMGh / k3 * (1 - Math.exp(-k3 * wearH));
                    return amtAtRemoval * Math.exp(-k3 * (tau - wearH));
                }
            }
            const paramsK = { ...params, k3 };
            const amtUnder = oneCompAmount(tau, event.doseMG, paramsK);
            if (tau > wearH) {
                const amtAtRemoval = oneCompAmount(wearH, event.doseMG, paramsK);
                return amtAtRemoval * Math.exp(-k3 * (tau - wearH));
            }
            return amtUnder;
        }
        default:
            return 0;
    }
}

/**
 * Compute CPA amount in central compartment (mg) for a single dose event,
 * scaling the effective dose by an adherence factor inferred from E2 EKF.
 *
 * Scientific basis: E2 θₛ (amplitude) encodes how faithfully actual medication
 * intake matches records.  This adherence factor is applied to the CPA dose,
 * NOT to CPA PK parameters — because E2 and CPA share the same dosing schedule
 * but their clearance mechanisms are independent (CPA is not directly measurable).
 *
 * @param adherenceScale  exp(θₛ) from E2 EKF (dose multiplier, not PK parameter)
 */
function computeCPAAmountWithAdherence(
    event: DoseEvent,
    tau: number,
    adherenceScale: number
): number {
    if (tau < 0 || event.ester !== Ester.CPA) return 0;
    return compute2CompCPACentralAmount(event.doseMG * adherenceScale, tau);
}

/**
 * Compute CPA plasma concentration (ng/mL) at a single time point using the
 * 2-compartment model, with E2-inferred adherence applied to the CPA dose.
 *
 * theta[0] (θₛ, amplitude) → adherenceScale = exp(θₛ): scales CPA effective dose.
 *   E2 amplitude encodes how closely recorded dosing matches actual intake.
 *   Applying this to CPA dose (not CPA PK) is scientifically defensible because:
 *   – E2 and CPA share the same dosing schedule/adherence behaviour.
 *   – CPA PK parameters (CL, Vd) are independent of E2 and cannot be inferred
 *     from E2 measurements alone.
 *
 * theta[1] (θₖ, clearance) is NOT applied to CPA — there are no CPA measurements
 * to constrain CPA-specific clearance variation.
 */
export function computeCPAAtTimeWithTheta(
    events: DoseEvent[],
    weight: number,
    timeH: number,
    theta: [number, number]
): number {
    const adherence = Math.exp(theta[0]); // dose-level adherence proxy from E2 EKF
    const sorted = [...events].sort((a, b) => a.timeH - b.timeH);
    let totalCentralMG = 0;
    for (const ev of sorted) {
        if (ev.timeH > timeH) continue;
        totalCentralMG += computeCPAAmountWithAdherence(ev, timeH - ev.timeH, adherence);
    }
    const V1_mL = CPA_2COMP_PK.V1_per_kg * weight * 1000;
    return Math.max(0, (totalCentralMG * 1e6) / V1_mL); // ng/mL
}

/**
 * Compute E2 plasma concentration (pg/mL) at a single time point,
 * applying individual parameters theta = [theta_s, theta_k].
 *
 * C(t; θ) = exp(θ_s) × C_pk(t, exp(θ_k))
 */
export function computeE2AtTimeWithTheta(
    events: DoseEvent[],
    weight: number,
    timeH: number,
    theta: [number, number]
): number {
    const s = Math.exp(theta[0]);
    const kScale = Math.exp(theta[1]);
    const sorted = [...events].sort((a, b) => a.timeH - b.timeH);

    let totalMG = 0;
    for (const ev of sorted) {
        if (ev.timeH > timeH) continue;
        totalMG += computeEventAmountWithKScale(ev, sorted, timeH - ev.timeH, kScale);
    }

    const plasmaVolML = CorePK.vdPerKG * weight * 1000;
    return Math.max(0, (totalMG * 1e9) / plasmaVolML * s);
}

/**
 * EKF update: incorporate one new lab result into the personal model.
 * Returns the updated state and diagnostics for that observation.
 */
export function ekfUpdatePersonalModel(
    events: DoseEvent[],
    weight: number,
    state: PersonalModelState,
    labResult: LabResult,
    prevLabTimeH?: number   // previous lab's timeH for time-scaled Q drift
): { newState: PersonalModelState; diagnostics: EKFDiagnostics } {
    const hasDoseBeforeLab = events.some(ev =>
        ev.timeH <= labResult.timeH &&
        ev.route !== Route.patchRemove &&
        ev.ester !== Ester.CPA
    );

    // --- Observation ---
    const obsPGmL = convertToPgMl(labResult.concValue, labResult.unit);
    const y = Math.log(Math.max(obsPGmL, EKF_EPS));

    // --- Prediction step (time-scaled parameter random walk) ---
    // Q is scaled by elapsed time relative to EKF_Q_REF_PERIOD_H (30 days).
    // This prevents over-shrinkage on long gaps and over-inflation on short gaps.
    const theta = state.thetaMean.slice() as [number, number];
    const dtH = (prevLabTimeH !== undefined)
        ? Math.max(24, labResult.timeH - prevLabTimeH)
        : EKF_Q_REF_PERIOD_H;
    const qScale = dtH / EKF_Q_REF_PERIOD_H;
    const P: [[number, number], [number, number]] = [
        [state.thetaCov[0][0] + state.Q[0][0] * qScale, state.thetaCov[0][1] + state.Q[0][1] * qScale],
        [state.thetaCov[1][0] + state.Q[1][0] * qScale, state.thetaCov[1][1] + state.Q[1][1] * qScale],
    ];

    // --- Predicted observation ---
    const predPGmL = computeE2AtTimeWithTheta(events, weight, labResult.timeH, theta);
    const yhat = Math.log(Math.max(predPGmL, EKF_EPS));

    // If there is no dosing history before this lab, treat it as a baseline point:
    // keep parameters unchanged and avoid flagging it as an outlier.
    if (!hasDoseBeforeLab) {
        const initialTrace = EKF_INITIAL_COV[0][0] + EKF_INITIAL_COV[1][1];
        const currentTrace = state.thetaCov[0][0] + state.thetaCov[1][1];
        const convergenceScore = Math.max(0, Math.min(1, 1 - currentTrace / initialTrace));

        const baselineState: PersonalModelState = {
            ...state,
            // Baseline points are recorded but do not contribute to EKF learning count.
            observationCount: state.observationCount,
            updatedAt: new Date().toISOString(),
        };

        const diagnostics: EKFDiagnostics = {
            NIS: 0,
            isOutlier: false,
            residualLog: 0,
            predictedPGmL: predPGmL,
            observedPGmL: obsPGmL,
            ci95Low: obsPGmL,
            ci95High: obsPGmL,
            convergenceScore,
            thetaS: Math.exp(state.thetaMean[0]),
            thetaK: Math.exp(state.thetaMean[1]),
        };

        return { newState: baselineState, diagnostics };
    }

    // --- Jacobian H = [∂yhat/∂theta_s, ∂yhat/∂theta_k] ---
    // ∂yhat/∂theta_s = 1 (exact, since yhat = theta_s + log(C_pk))
    const thetaKPerturbed: [number, number] = [theta[0], theta[1] + EKF_DELTA_K];
    const predPerturbed = computeE2AtTimeWithTheta(events, weight, labResult.timeH, thetaKPerturbed);
    const yhatPerturbed = Math.log(Math.max(predPerturbed, EKF_EPS));
    const H: [number, number] = [1.0, (yhatPerturbed - yhat) / EKF_DELTA_K];

    // --- Innovation ---
    const nu = y - yhat;

    // --- Innovation covariance S = H P H^T + R ---
    const S = H[0]*H[0]*P[0][0] + 2*H[0]*H[1]*P[0][1] + H[1]*H[1]*P[1][1] + state.Rlog;

    // --- Outlier detection (Normalised Innovation Squared) ---
    const NIS = (S > 0) ? (nu * nu / S) : 0;
    const isOutlier = NIS > EKF_CHI2_95;

    // --- Inflate R if outlier (robust update) ---
    const Reff = isOutlier ? state.Rlog * 4.0 : state.Rlog;
    const Seff = H[0]*H[0]*P[0][0] + 2*H[0]*H[1]*P[0][1] + H[1]*H[1]*P[1][1] + Reff;

    // --- Kalman gain K = P H^T / Seff ---
    const K: [number, number] = [
        (P[0][0]*H[0] + P[0][1]*H[1]) / Seff,
        (P[1][0]*H[0] + P[1][1]*H[1]) / Seff,
    ];

    // --- Update theta ---
    const thetaNew: [number, number] = [theta[0] + K[0] * nu, theta[1] + K[1] * nu];

    // --- Update covariance  P_new = (I - K H) P ---
    const i00 = 1 - K[0]*H[0];
    const i01 = -K[0]*H[1];
    const i10 = -K[1]*H[0];
    const i11 = 1 - K[1]*H[1];
    const PNew: [[number, number], [number, number]] = [
        [i00*P[0][0] + i01*P[1][0], i00*P[0][1] + i01*P[1][1]],
        [i10*P[0][0] + i11*P[1][0], i10*P[0][1] + i11*P[1][1]],
    ];
    // Enforce symmetry and positive lower bounds
    PNew[0][1] = PNew[1][0] = (PNew[0][1] + PNew[1][0]) / 2;
    PNew[0][0] = Math.max(PNew[0][0], 1e-6);
    PNew[1][1] = Math.max(PNew[1][1], 1e-6);

    // --- Residual anchor (log-space residual after update) ---
    const newPredPGmL = computeE2AtTimeWithTheta(events, weight, labResult.timeH, thetaNew);
    const logRatioPost = Math.log(Math.max(obsPGmL, EKF_EPS)) - Math.log(Math.max(newPredPGmL, EKF_EPS));
    const anchor: ResidualAnchor = {
        timeH: labResult.timeH,
        logRatio: logRatioPost,
        w: isOutlier ? 0.3 : 1.0,
        kind: 'lab',
    };
    const updatedAnchors = [...state.anchors, anchor]
        .sort((a, b) => a.timeH - b.timeH)
        .slice(-20); // keep most recent 20

    // --- 95% CI at observation time ---
    const hK = H[1];
    const varYhat = PNew[0][0] + 2*PNew[0][1]*hK + PNew[1][1]*hK*hK;
    const std95 = Math.sqrt(Math.max(0, varYhat + Reff));
    const logPredNew = Math.log(Math.max(newPredPGmL, EKF_EPS));
    const ci95Low = Math.exp(logPredNew - 1.96 * std95);
    const ci95High = Math.exp(logPredNew + 1.96 * std95);

    // --- Convergence score (how much uncertainty has been reduced) ---
    const initialTrace = EKF_INITIAL_COV[0][0] + EKF_INITIAL_COV[1][1];
    const currentTrace = PNew[0][0] + PNew[1][1];
    const convergenceScore = Math.max(0, Math.min(1, 1 - currentTrace / initialTrace));

    const newState: PersonalModelState = {
        modelVersion: 'pk-ekf-v1',
        thetaMean: thetaNew,
        thetaCov: PNew,
        Q: state.Q,
        Rlog: state.Rlog,
        anchors: updatedAnchors,
        observationCount: state.observationCount + 1,
        updatedAt: new Date().toISOString(),
    };

    const diagnostics: EKFDiagnostics = {
        NIS,
        isOutlier,
        residualLog: nu,
        predictedPGmL: predPGmL,
        observedPGmL: obsPGmL,
        ci95Low,
        ci95High,
        convergenceScore,
        thetaS: Math.exp(thetaNew[0]),
        thetaK: Math.exp(thetaNew[1]),
    };

    return { newState, diagnostics };
}

/**
 * Replay all lab results from prior state to rebuild the personal model.
 * Call this when events are edited/deleted or lab results are changed.
 */
export function replayPersonalModel(
    events: DoseEvent[],
    weight: number,
    labResults: LabResult[]
): PersonalModelState {
    let state = initPersonalModel();
    const sorted = [...labResults].sort((a, b) => a.timeH - b.timeH);
    for (let i = 0; i < sorted.length; i++) {
        const prevTimeH = i > 0 ? sorted[i - 1].timeH : undefined;
        const { newState } = ekfUpdatePersonalModel(events, weight, state, sorted[i], prevTimeH);
        state = newState;
    }
    return state;
}

/**
 * Compute the fraction of E2 clearance inhibited by CPA at a given concentration.
 *
 * Uses the Chou–Talalay / Hill median-effect model:
 *   f_a = Emax · CPA^n / (IC50^n + CPA^n)
 *
 * Returns a value in [0, Emax]. Pass into computeSimulationWithCI via
 * applyCPAInhibitionToE2 = true to apply the post-hoc E2 scaling.
 *
 * @param cpaNgML  CPA plasma concentration in ng/mL
 * @param params   Hill model parameters (defaults to CPA_E2_INHIBITION)
 */
export function computeCPAE2InhibitionFactor(
    cpaNgML: number,
    params = CPA_E2_INHIBITION
): number {
    if (!Number.isFinite(cpaNgML) || cpaNgML <= 0) return 0;
    const { Emax, IC50, n } = params;
    if (!Number.isFinite(Emax) || Emax <= 0 || !Number.isFinite(IC50) || IC50 <= 0 || !Number.isFinite(n) || n <= 0) return 0;
    const Dn    = Math.pow(cpaNgML, n);
    const IC50n = Math.pow(IC50, n);
    // Clamp result to [0, Emax) to guarantee scale = 1/(1-f) is always well-defined
    return Math.min(Math.max(0, Emax * Dn / (IC50n + Dn)), Emax * 0.9999);
}

/**
 * Compute a full simulation curve with Bayesian OU-Kalman calibration for E2
 * and EKF-based adherence scaling for CPA.
 *
 * E2 calibration uses either:
 *  - 'ekf': EKF-learned (theta_s, theta_k) parameters applied to the PK model;
 *           CI from linearised parameter uncertainty H^T P H.
 *  - 'ou-kalman': OU Kalman filter + RTS smoother log-factor x(t);
 *                 c(t) = c₀(t)·exp(x(t)), bands at ±σ (68%) and ±1.96σ (95%).
 *
 * CPA adherence scaling uses EKF theta[0] in both modes.
 *
 * @param labResults            Raw lab results needed by the OU-Kalman filter.
 * @param calibrationModel      Which E2 model to use (default 'ekf').
 * @param applyCPAInhibitionToE2 When true, applies Chou-Talalay CPA→E2 clearance
 *                               inhibition as a post-hoc multiplicative correction
 *                               (default false — experimental, off until validated).
 */
export function computeSimulationWithCI(
    sim: SimulationResult,
    events: DoseEvent[],
    weight: number,
    state: PersonalModelState,
    applyE2LearningToCPA: boolean = true,
    labResults: LabResult[] = [],
    calibrationModel: CalibrationModel = 'ekf',
    applyCPAInhibitionToE2: boolean = false
): {
    timeH: number[];
    e2Adjusted: number[];
    ci95Low: number[];
    ci95High: number[];
    ci68Low: number[];
    ci68High: number[];
    cpaAdjusted: number[];
    cpaCi95Low: number[];
    cpaCi95High: number[];
} {
    const n = sim.timeH.length;
    if (n === 0) {
        return { timeH: [], e2Adjusted: [], ci95Low: [], ci95High: [], ci68Low: [], ci68High: [], cpaAdjusted: [], cpaCi95Low: [], cpaCi95High: [] };
    }
    const theta = state.thetaMean;
    const P = state.thetaCov;

    const clampCI = (low: number, high: number, hardMax: number): [number, number] => {
        const lo = Number.isFinite(low) ? Math.max(0, low) : 0;
        const hi = Number.isFinite(high) ? Math.max(lo, high) : lo;
        return [Math.min(lo, hardMax), Math.min(hi, hardMax)];
    };

    // ── E2: choose calibration model ──────────────────────────────────────────
    const e2Adjusted = new Array<number>(n);
    const ci95Low    = new Array<number>(n);
    const ci95High   = new Array<number>(n);
    const ci68Low    = new Array<number>(n);
    const ci68High   = new Array<number>(n);

    if (calibrationModel === 'ou-kalman') {
        // OU-Kalman: posterior mean m(t) and variance P(t) from forward Kalman + RTS smoother.
        // c(t) = c₀(t) · exp(m(t))  [median prediction]
        // CI bands are symmetric in log-space → lognormal in concentration space.
        const ou = buildOUKalmanCalibration(sim, labResults);

        for (let i = 0; i < n; i++) {
            const c0  = Math.max(sim.concPGmL_E2[i], EKF_EPS);
            const mi  = ou.m[i];
            const std = Math.sqrt(Math.max(0, ou.P[i]));

            // Jensen mean correction: E[c] = c₀·exp(m + P/2) (lognormal mean > median)
            e2Adjusted[i] = Math.min(c0 * Math.exp(mi + 0.5 * ou.P[i]), EKF_CI_MAX_E2);

            const [lo95, hi95] = clampCI(
                c0 * Math.exp(mi - 1.96 * std),
                c0 * Math.exp(mi + 1.96 * std),
                EKF_CI_MAX_E2
            );
            ci95Low[i]  = lo95;
            ci95High[i] = hi95;

            const [lo68, hi68] = clampCI(
                c0 * Math.exp(mi - std),
                c0 * Math.exp(mi + std),
                EKF_CI_MAX_E2
            );
            ci68Low[i]  = lo68;
            ci68High[i] = hi68;
        }
    } else {
        // EKF: use learned (theta_s, theta_k) parameters to shift the PK curve.
        // CI from linearised parameter uncertainty H^T P H (where H=[1, ∂logC/∂theta_k]).
        // Central curve applies Jensen mean correction: E[C] = C_MAP·exp(σ²_param/2),
        // which compensates for the lognormal mean > median bias under parameter uncertainty.
        // Sampled at ~100 points and linearly interpolated.
        const ekfStep = Math.max(1, Math.floor(n / 100));
        const ekfIndices: number[] = [];
        for (let i = 0; i < n; i += ekfStep) ekfIndices.push(i);
        if (ekfIndices[ekfIndices.length - 1] !== n - 1) ekfIndices.push(n - 1);

        interface EKFSample {
            idx: number;
            e2Adj: number;   // Jensen mean-corrected central curve (E[C])
            ci95Lo: number;  // 5th percentile of posterior predictive (anchored to MAP)
            ci95Hi: number;  // 95th percentile
            ci68Lo: number;  // 16th percentile
            ci68Hi: number;  // 84th percentile
        }
        const ekfSamples: EKFSample[] = [];

        for (const idx of ekfIndices) {
            const timeH  = sim.timeH[idx];
            const e2Base = computeE2AtTimeWithTheta(events, weight, timeH, theta);
            const yhat   = Math.log(Math.max(e2Base, EKF_EPS));
            // Jacobian H[1] = ∂log(C)/∂theta_k via forward finite difference
            const thetaKPlus: [number, number] = [theta[0], theta[1] + EKF_DELTA_K];
            const yhatPlus = Math.log(Math.max(computeE2AtTimeWithTheta(events, weight, timeH, thetaKPlus), EKF_EPS));
            const H1 = (yhatPlus - yhat) / EKF_DELTA_K;
            // σ²_param = H^T P H (parameter uncertainty only)
            const rawSigma2Param = P[0][0] + 2 * H1 * P[0][1] + H1 * H1 * P[1][1];
            const sigma2Param = (Number.isFinite(rawSigma2Param) && rawSigma2Param > 0) ? rawSigma2Param : 0;
            // σ²_total = σ²_param + σ²_residual (add irreducible assay + biological noise)
            const sigmaTotal = Math.sqrt(sigma2Param + EKF_SIGMA_RESIDUAL_LOG * EKF_SIGMA_RESIDUAL_LOG);
            // Jensen mean correction: E[C] = C_MAP × exp(σ²_param / 2)
            const e2AdjMean = Math.min(e2Base * Math.exp(0.5 * sigma2Param), EKF_CI_MAX_E2);
            // CI percentile bounds are anchored to the MAP (lognormal median)
            const [lo95, hi95] = clampCI(e2Base * Math.exp(-1.96 * sigmaTotal), e2Base * Math.exp(1.96 * sigmaTotal), EKF_CI_MAX_E2);
            const [lo68, hi68] = clampCI(e2Base * Math.exp(-sigmaTotal), e2Base * Math.exp(sigmaTotal), EKF_CI_MAX_E2);
            ekfSamples.push({ idx, e2Adj: e2AdjMean, ci95Lo: lo95, ci95Hi: hi95, ci68Lo: lo68, ci68Hi: hi68 });
        }

        for (let j = 0; j < ekfSamples.length; j++) {
            const a = ekfSamples[j];
            const b = ekfSamples[j + 1] ?? a;
            const span = b.idx - a.idx;
            for (let i = a.idx; i <= b.idx; i++) {
                const frac = span > 0 ? (i - a.idx) / span : 0;
                e2Adjusted[i] = Math.min(Math.max(0, a.e2Adj + (b.e2Adj - a.e2Adj) * frac), EKF_CI_MAX_E2);
                ci95Low[i]  = a.ci95Lo + (b.ci95Lo - a.ci95Lo) * frac;
                ci95High[i] = a.ci95Hi + (b.ci95Hi - a.ci95Hi) * frac;
                ci68Low[i]  = a.ci68Lo + (b.ci68Lo - a.ci68Lo) * frac;
                ci68High[i] = a.ci68Hi + (b.ci68Hi - a.ci68Hi) * frac;
            }
        }
    }

    // ── CPA: EKF adherence-scaled prediction + population PK uncertainty ──────
    // theta[0] (E2 amplitude) proxies adherence for CPA dose scaling.
    // theta[1] (clearance) is NOT applied to CPA — no CPA blood measurements.
    const step = Math.max(1, Math.floor(n / 100));
    const sampledIndices: number[] = [];
    for (let i = 0; i < n; i += step) sampledIndices.push(i);
    if (sampledIndices[sampledIndices.length - 1] !== n - 1) sampledIndices.push(n - 1);

    const cpaSamples: {
        idx: number;
        cpaAdj: number;
        cpaCi95Low: number;
        cpaCi95High: number;
    }[] = [];

    for (const idx of sampledIndices) {
        const timeH = sim.timeH[idx];

        const cpaPred = applyE2LearningToCPA
            ? computeCPAAtTimeWithTheta(events, weight, timeH, theta)   // adherence-scaled
            : computeCPAAtTimeWithTheta(events, weight, timeH, [0, 0]); // population mean

        // CPA CI: adherence uncertainty (P[0][0] when enabled) + population PK variance.
        const adherenceVar = applyE2LearningToCPA ? Math.max(0, P[0][0]) : 0;
        const varLogCPA    = adherenceVar + CPA_2COMP_PK.popLogVar;
        const stdCPA       = Math.sqrt(Math.max(0, varLogCPA));
        const yhatCPA      = Math.log(Math.max(cpaPred, EKF_EPS_CPA));
        const [cpaCiLow, cpaCiHigh] = clampCI(
            Math.exp(yhatCPA - 1.96 * stdCPA),
            Math.exp(yhatCPA + 1.96 * stdCPA),
            EKF_CI_MAX_CPA
        );

        cpaSamples.push({ idx, cpaAdj: cpaPred, cpaCi95Low: cpaCiLow, cpaCi95High: cpaCiHigh });
    }

    // Linear interpolation for CPA across all n points
    const cpaAdjusted = new Array<number>(n).fill(0);
    const cpaCi95Low  = new Array<number>(n).fill(0);
    const cpaCi95High = new Array<number>(n).fill(0);

    for (let j = 0; j < cpaSamples.length; j++) {
        const a = cpaSamples[j];
        const b = cpaSamples[j + 1] ?? a;
        const span = b.idx - a.idx;
        for (let i = a.idx; i <= b.idx; i++) {
            const frac = span > 0 ? (i - a.idx) / span : 0;
            cpaAdjusted[i]  = a.cpaAdj     + (b.cpaAdj     - a.cpaAdj)     * frac;
            cpaCi95Low[i]   = a.cpaCi95Low  + (b.cpaCi95Low  - a.cpaCi95Low)  * frac;
            cpaCi95High[i]  = a.cpaCi95High + (b.cpaCi95High - a.cpaCi95High) * frac;
        }
    }

    // ── Optional: CPA → E2 clearance inhibition (Chou-Talalay Hill model) ──────
    // Applies a post-hoc multiplicative correction to E2 and its CI bands.
    // Steady-state approximation: lower clearance → concentration ∝ 1 / (1 − f).
    // Default off (experimental); enable when supported by user's CPA blood levels.
    if (applyCPAInhibitionToE2) {
        for (let i = 0; i < n; i++) {
            // Use adherence-adjusted CPA curve (consistent with displayed CPA trajectory)
            const f     = computeCPAE2InhibitionFactor(cpaAdjusted[i]);
            const scale = 1 / (1 - f); // f ∈ [0, 0.25] → scale ∈ [1.0, 1.333]
            e2Adjusted[i]  = Math.min(e2Adjusted[i]  * scale, EKF_CI_MAX_E2);
            // Scale all CI bounds then re-clamp both lower and upper to avoid ordering violations
            const s95Lo = ci95Low[i]  * scale;
            const s95Hi = Math.min(ci95High[i] * scale, EKF_CI_MAX_E2);
            const s68Lo = ci68Low[i]  * scale;
            const s68Hi = Math.min(ci68High[i] * scale, EKF_CI_MAX_E2);
            ci95Low[i]  = Math.min(s95Lo, s95Hi);
            ci95High[i] = s95Hi;
            ci68Low[i]  = Math.min(s68Lo, s68Hi);
            ci68High[i] = s68Hi;
        }
    }

    return {
        timeH: sim.timeH,
        e2Adjusted, ci95Low, ci95High, ci68Low, ci68High,
        cpaAdjusted, cpaCi95Low, cpaCi95High,
    };
}

async function generateKey(password: string, salt: Uint8Array) {
    const enc = new TextEncoder();
    const keyMaterial = await window.crypto.subtle.importKey(
        "raw",
        enc.encode(password),
        { name: "PBKDF2" },
        false,
        ["deriveKey"]
    );
    return window.crypto.subtle.deriveKey(
        {
            name: "PBKDF2",
            salt: salt as any,
            iterations: 100000,
            hash: "SHA-256"
        },
        keyMaterial,
        { name: "AES-GCM", length: 256 },
        false,
        ["encrypt", "decrypt"]
    );
}

function buffToBase64(buff: Uint8Array): string {
    const bin = Array.from(buff, (byte) => String.fromCharCode(byte)).join("");
    return btoa(bin);
}

function base64ToBuff(b64: string): Uint8Array {
    const bin = atob(b64);
    return Uint8Array.from(bin, (c) => c.charCodeAt(0));
}

export async function encryptData(text: string): Promise<{ data: string, password: string }> {
    const password = buffToBase64(window.crypto.getRandomValues(new Uint8Array(12)));
    const salt = window.crypto.getRandomValues(new Uint8Array(16));
    const iv = window.crypto.getRandomValues(new Uint8Array(12));
    const key = await generateKey(password, salt);
    const enc = new TextEncoder();
    const encrypted = await window.crypto.subtle.encrypt(
        { name: "AES-GCM", iv: iv as any },
        key,
        enc.encode(text)
    );

    const bundle = {
        encrypted: true,
        iv: buffToBase64(iv),
        salt: buffToBase64(salt),
        data: buffToBase64(new Uint8Array(encrypted))
    };
    return {
        data: JSON.stringify(bundle),
        password
    };
}

export async function decryptData(jsonString: string, password: string): Promise<string | null> {
    try {
        const bundle = JSON.parse(jsonString);
        if (!bundle.encrypted) return jsonString;

        const salt = base64ToBuff(bundle.salt);
        const iv = base64ToBuff(bundle.iv);
        const data = base64ToBuff(bundle.data);

        const key = await generateKey(password, salt);
        const decrypted = await window.crypto.subtle.decrypt(
            { name: "AES-GCM", iv: iv as any },
            key,
            data as any
        );
        return new TextDecoder().decode(decrypted);
    } catch (e) {
        console.error(e);
        return null;
    }
}
