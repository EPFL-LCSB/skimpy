def stringify_stoichiometry(stoichiometry):
    """
    Stoichiometry looks like [-2, -1, 3]
    This is not nice as a suffix for a class name
    We transform it to _m2_m1_3
    :param stoichiometry:
    :type stoichiometry: list(float)
    :return: suffix
    """

    suffix = '_'.join([str(x) for x in stoichiometry])
    suffix = suffix.replace('-','m')
    return suffix