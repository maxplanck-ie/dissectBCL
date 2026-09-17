import pandas as pd


def parseUnclassifiedPct(reportPath):
    """
    Returns the % of reads reported as 'unclassified' in a kraken2
    --report TSV, or None if the report is empty/unreadable, or has
    no unclassified row.
    """
    try:
        reportDF = pd.read_csv(reportPath, sep="\t", header=None)
        unclassified = reportDF[reportDF[5].str.strip() == "unclassified"]
    except (pd.errors.EmptyDataError, pd.errors.ParserError, KeyError, OSError):
        return None
    if unclassified.empty:
        return None
    return float(unclassified.iloc[0][0])


def pickThreshold(analysisType, config):
    """
    Returns the unclassified-% threshold (a float) that applies to
    analysisType: config['screening']['relaxed_threshold'] if analysisType
    is listed in config['screening']['relaxed_analysis_types'] (matched
    case-insensitively), else config['screening']['unclassified_threshold'].
    """
    relaxedTypes = {
        t.strip().lower()
        for t in config["screening"]
        .get("relaxed_analysis_types", fallback="")
        .split(",")
        if t.strip()
    }
    if isinstance(analysisType, str) and analysisType.strip().lower() in relaxedTypes:
        return config["screening"].getfloat("relaxed_threshold", fallback=10.0)
    return config["screening"].getfloat("unclassified_threshold", fallback=10.0)


def needsEscalation(reportPath, analysisType, config):
    """
    Returns True if reportPath's unclassified % exceeds the threshold
    that applies to analysisType (see pickThreshold).
    """
    pct = parseUnclassifiedPct(reportPath)
    if pct is None:
        return False
    return pct > pickThreshold(analysisType, config)
