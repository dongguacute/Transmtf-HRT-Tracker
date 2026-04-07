// Re-export the shared domain model from dedicated files so existing imports
// from `logic.ts` keep working while we progressively untangle responsibilities.
export { Route, Ester, ExtraKey, type DoseEvent, type SimulationResult, type LabResult } from './types';
export {
    getToE2Factor,
    SL_TIER_ORDER,
    SublingualTierParams,
    getBioavailabilityMultiplier,
    runSimulation,
    interpolateConcentration,
    interpolateConcentration_E2,
    interpolateConcentration_CPA,
} from './pk';
export {
    convertToPgMl,
    createCalibrationInterpolator,
    OU_DEFAULT_PARAMS,
    buildOUKalmanCalibration,
    type CalibrationModel,
    type OUCalibParams,
} from './calibration';
export {
    initPersonalModel,
    computeCPAAtTimeWithTheta,
    computeE2AtTimeWithTheta,
    ekfUpdatePersonalModel,
    replayPersonalModel,
    computeCPAE2InhibitionFactor,
    computeSimulationWithCI,
    type ResidualAnchor,
    type PersonalModelState,
    type EKFDiagnostics,
} from './personalModel';
export { encryptData, decryptData } from './src/utils/dataEncryption';
