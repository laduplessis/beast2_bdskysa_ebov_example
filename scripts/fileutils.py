import yaml





def parseList(clist, datatype=str, sep=","):
    """Parse a comma-separated input argument as a list of the given datatype

    datatype is a function that is applied to the list entries

    """

    if (clist == ""):
    	return []
    else:
	    return list(map(datatype,clist.split(sep)))




def getMapping(header, fieldmapping):
    """Figure out the mapping between the database columns and columns we are processing

    """ 
    
    mapping = dict()
    for field in fieldmapping.keys():        
        if (isinstance(fieldmapping[field], str)):
                fieldmapping[field] = [fieldmapping[field]]

        for i in range(0,len(header)):                    
            if (header[i].strip() in fieldmapping[field]):
                mapping[field] = i

    return mapping
#

