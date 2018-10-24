import sys
import time

def read_seq_line(name):
    fp = open(name)
    lines = fp.readlines()
    fp.close()

    return lines[1]

def read_score_file(name):
    fp = open(name)
    lines = fp.readlines()
    fp.close()
    
    match = int(lines[0].split('=')[1])
    mismatch = int(lines[1].split('=')[1])
    gap = int(lines[2].split('=')[1])
    
    return match, mismatch, gap
        
def global_allign(seq1, seq2):
    '''
      seq1      : sequence1 string value 
      seq2      : sequence2 string value

      mat       : scores will be stored
      backtrack : tracks will be stored
    '''

    mat = [[None for x in range(len(seq2))] for y in range(len(seq1))]
    backtrack =  [[None for x in range(len(seq2))] for y in range(len(seq1))]

    mat[0][0] = backtrack[0][0] = 0
    
    for i in range(len(seq1)):
        for j in range(len(seq2)):
            if(mat[i][j] is None):
                if i == 0:
                    mat[i][j] = mat[i][j-1] + gap
                    backtrack[i][j] = 'left'
                elif j == 0:
                    mat[i][j] = mat[i-1][j] + gap
                    backtrack[i][j] = 'up'
                else:
                    mat[i][j] = max(mat[i][j-1] + gap, 
                                    mat[i-1][j] + gap,
                                    mat[i-1][j-1] + (match if (seq1[i] == seq2[j]) else mismatch))
                    
                    if mat[i][j] == mat[i][j-1] + gap:
                        backtrack[i][j] = 'left'
                    elif mat[i][j] == mat[i-1][j] + gap:
                        backtrack[i][j] = 'up'
                    else:
                        backtrack[i][j] = 'diag'
    
    # backtracking results will be stored in 'str1', 'str2' 
    str1 = ''
    str2 = ''

    i = len(seq1)-1
    j = len(seq2)-1
    while(1):
        # Loop from end point to start point
        if i == 0 and j == 0:
            return str1, str2, mat[len(seq1)-1][len(seq2)-1]
        elif backtrack[i][j] == 'left':
            str1 = '-' + str1
            str2 = seq2[j] + str2
            j = j-1
        elif backtrack[i][j] == 'up':
            str1 = seq1[i] + str1
            str2 = '-' + str2
            i = i-1
        elif backtrack[i][j] == 'diag':
            str1 = seq1[i] + str1
            str2 = seq2[j] + str2
            i = i-1
            j = j-1
        else:
            print('Something went wrong')

def local_allign(seq1, seq2):
    local_mat = [[None for x in range(len(seq2))] for y in range(len(seq1))]
    backtrack =  [[None for x in range(len(seq2))] for y in range(len(seq1))]
    
    local_max = 0
    max_i = max_j = 0
    
    local_mat[0][0] = backtrack[0][0] = 0
    
    for i in range(len(seq1)):
        for j in range(len(seq2)):
            if(local_mat[i][j] is None):
                if i == 0:
                    local_mat[i][j] = 0
                    backtrack[i][j] = 'left'
                elif j == 0:
                    local_mat[i][j] = 0
                    backtrack[i][j] = 'up'
                else:
                    local_mat[i][j] = max(0,
                                          local_mat[i][j-1] + gap, 
                                          local_mat[i-1][j] + gap,
                                          local_mat[i-1][j-1] + (match if (seq1[i] == seq2[j]) else mismatch))
                    
                    if local_mat[i][j] == local_mat[i][j-1] + gap:
                        backtrack[i][j] = 'left'
                    elif local_mat[i][j] == local_mat[i-1][j] + gap:
                        backtrack[i][j] = 'up'
                    else:
                        backtrack[i][j] = 'diag'
                        
                if(local_mat[i][j] > local_max):
                    local_max = local_mat[i][j]
                    max_i = i
                    max_j = j
    
    #BackTracking
    str1 = ''
    str2 = ''
   
    i = max_i
    j = max_j
    
    while(1):
        # Start from maximum score
        # Loop until reaches the point the score is zero
        if local_mat[i][j] == 0:
            return str1, str2, local_mat[max_i][max_j]
        elif backtrack[i][j] == 'left':
            str1 = '-' + str1
            str2 = seq2[j] + str2
            j = j-1
        elif backtrack[i][j] == 'up':
            str1 = seq1[i] + str1
            str2 = '-' + str2
            i = i-1
        elif backtrack[i][j] == 'diag':
            str1 = seq1[i] + str1
            str2 = seq2[j] + str2
            i = i-1
            j = j-1
        else:
            print('Something went wrong')
    
    return str1, str2, mat[len(seq1)-1][len(seq2)-1]


if __name__ == "__main__":
    file1 = sys.argv[1] + '.txt'
    file2 = sys.argv[2] + '.txt'
    scoreFile = sys.argv[3] + '.txt'

    # file1 = 'seq1.txt'
    # file2 = 'seq2.txt'
    # scoreFile = 'score.txt'

    # append space at the beginning
    seq1 = ' ' + read_seq_line(file1)
    seq2 = ' ' + read_seq_line(file2)
    match, mismatch, gap = read_score_file(scoreFile)

    g_str1, g_str2, g_score = global_allign(seq1, seq2)
    l_str1, l_str2, l_score = local_allign(seq1, seq2)

    print('Global Score\t: %d'%(g_score))
    print('Local Score\t: %d\n'%(l_score))

    cur_time = time.strftime("%y-%m-%d at %H.%M.%S")
    global_result_name =  'result_g ' + cur_time + '.txt'
    local_result_name = 'result_l ' + cur_time + '.txt'

    fp = open(global_result_name,'w')
    fp.write('Global_Seq1 : %s\nGlobal Seq2 : %s\nGlobal Score: %d'%(g_str1, g_str2, g_score))
    fp.close()
    
    fp = open(local_result_name,'w')
    fp.write('Local Seq1 : %s\nLocal Seq2 : %s\nLocal Score: %d'%(l_str1, l_str2, l_score))
    fp.close()
