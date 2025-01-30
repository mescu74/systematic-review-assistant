from __future__ import annotations

from enum import StrEnum


class PopulationExclusion(StrEnum):
    AGE_RANGE = "Wrong age range for inclusion criteria"
    SPECIES = "Non-human subjects"
    CONDITION = "Condition or disease mismatch"
    COMORBIDITY = "Excluded comorbidity present"
    DEMOGRAPHIC = "Demographics outside criteria"
    SETTING = "Inappropriate study setting"

class InterventionExclusion(StrEnum):
    TIMING = "Intervention timing does not match criteria"
    DOSAGE = "Dosage outside specified range"
    DELIVERY = "Incorrect delivery or administration method"
    PROTOCOL = "Protocol deviation from inclusion criteria"
    COMBINATION = "Excluded intervention combination"
    DURATION = "Intervention duration outside criteria"

class ComparisonExclusion(StrEnum):
    CONTROL_TYPE = "Inappropriate control group"
    BASELINE = "Missing baseline data"
    COMPARATOR = "Incorrect comparator"
    RANDOMIZATION = "Inadequate randomization method"
    BLINDING = "Insufficient blinding procedure"
    CROSSOVER = "Inappropriate crossover design"

class OutcomeExclusion(StrEnum):
    MEASUREMENT = "Invalid outcome measurement method"
    TIMEPOINT = "Incorrect measurement timepoint"
    METRIC = "Inappropriate outcome metric"
    REPORTING = "Incomplete outcome reporting"
    SURROGATE = "Invalid surrogate endpoint"
    FOLLOWUP = "Insufficient follow-up period"

class StudyDesignExclusion(StrEnum):
    STUDY_TYPE = "Wrong study design type"
    SAMPLE_SIZE = "Sample size below requirement"
    DURATION = "Study duration too short"
    ETHICS = "Ethical concerns identified"
    PROTOCOL = "Major protocol violations"
    REGISTRATION = "Study not pre-registered"

class ReportingExclusion(StrEnum):
    LANGUAGE = "Non-included language"
    DUPLICATE = "Duplicate publication"
    DATE_RANGE = "Outside date range"
    PEER_REVIEW = "Not peer-reviewed"
    #ABSTRACT_ONLY = "Abstract only available"
    RETRACTION = "Retracted publication"

type ExclusionReasonType = (
    PopulationExclusion |
    InterventionExclusion |
    ComparisonExclusion |
    OutcomeExclusion |
    StudyDesignExclusion |
    ReportingExclusion
)
