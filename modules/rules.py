def all_rules(infoAln, dend):
	"""Set of all rules"""
	return any([ref_inside(infoAln, dend), read_inside(infoAln, dend), read_left_aligned(infoAln, dend), read_right_aligned(infoAln, dend)])


def ref_inside(infoAln, dend):
    """ 1. alignment over
	 all read length and ref length
	################################
	# REF     **********
	# aln     |||||||||| 
	# read *****************
	################################
	"""
    if (infoAln.targetStart <= dend and (infoAln.targetLength - dend) <= infoAln.targetEnd):
        return True


def read_inside(infoAln, dend):
    """ 2. alignment over all read length but on a small fraction of the ref
	################################
	# REF  **********
	# aln     |||| 
	# read    ****
	################################
	"""
    if (infoAln.queryStart <= dend and (infoAln.queryLength - dend) <= infoAln.queryEnd):
        return True


def read_right_aligned(infoAln, dend):
    """ 3. alignment half the read to the read right end
	################################
	# REF  		  **********
	# aln         |||||| 
	# read *************
	################################
	"""
    if (infoAln.targetStart <= dend and (infoAln.queryLength - dend) <= infoAln.queryEnd):
        return True


def read_left_aligned(infoAln, dend):
    """ 4. alignment from the read left end to the half of the read
	################################
	# REF  **********
	# aln      |||||| 
	# read     *************
	################################
	"""
    if (infoAln.queryStart <= dend and (infoAln.targetLength - dend) <= infoAln.targetEnd):
        return True
