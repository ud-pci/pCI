

def calc_lifetimes(tr_file):
    
    # Read transition rates from file tr_file
    f = open(tr_file, 'r')
    lines = f.readlines()
    f.close()
    
    # Parse transition rates file for energy levels and sum up transition rates
    tr_rates = {}
    for line in lines[1:]:
        configuration1 = line.split(',')[0]
        term1 = line.split(',')[1]
        J1 = line.split(',')[2]
        configuration2 = line.split(',')[3]
        term2 = line.split(',')[4]
        J2 = line.split(',')[5]
        tr_rate = float(line.split(',')[8][:-1])
        try:
            tr_rates[configuration2].append(tr_rate)
        except KeyError:
            tr_rates[configuration2] = []
            tr_rates[configuration2].append(tr_rate)
    
    for config, rates in tr_rates.items():
        lifetime = 1/sum(rates)
        print(config, lifetime)

if __name__ == "__main__":
    calc_lifetimes('tr_test.csv')