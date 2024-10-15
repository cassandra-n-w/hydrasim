def load_data_as_dic(filepath, comments='!', returnOriginalKeys=False):
    """load iter_000*.dat into a dictionary 
       from Fujun Du's pylib"""

    data = np.loadtxt(filepath, comments=comments)

    ftmp = open(filepath, 'r')
    str_comment = ftmp.readline()[1:].split()
    ftmp.close()

    dic = {}
    for i in xrange(len(str_comment)):
        dic.update({str_comment[i]: data[:, i]})

    del data

    if returnOriginalKeys:
        return str_comment, dic
    else:
        return dic
