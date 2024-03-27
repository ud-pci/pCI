
def set_min_uncertainties(filename, min_percentage):
    
    f = open(filename,'r')
    lines = f.readlines()
    f.close()
    
    matrix_res = []
    for line in lines:
        matrix_res.append(line.replace('\n','').split(','))
    
    cnt = 0
    for line in matrix_res[1:]:
        val = float(line[6]) # matrix element value
        try:
            unc = float(line[7]) # matrix element uncertainty
        except ValueError:
            continue
        if val*min_percentage/100 >= unc: cnt += 1
        line[7] = '{:,.5f}'.format(max(val*min_percentage/100, unc))
        
    new_filename = filename.split('.')[0] + '_adjusted.' + filename.split('.')[1]
    print('# of matrix element uncertainties adjusted:', cnt)
    f = open(new_filename, 'w')
    for line in matrix_res:
        f.write(','.join(line) + '\n')
    f.close()
    print(new_filename + ' has been written')

if __name__ == '__main__':
    filename = input('Input name of matrix element csv file: ')
    min_percent = float(input('Input minimum uncertainty to set in percentage: '))
    set_min_uncertainties(filename, min_percent)
    
    