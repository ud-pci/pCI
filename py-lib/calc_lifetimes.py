

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
        termJ1 = term1 + J1
        configuration2 = line.split(',')[3]
        term2 = line.split(',')[4]
        J2 = line.split(',')[5]
        termJ2 = term2 + J2
        tr_rate = float(line.split(',')[8][:-1])
        config2 = configuration2 + ' ' + termJ2
        try:
            tr_rates[config2].append(tr_rate)
        except KeyError:
            tr_rates[config2] = []
            tr_rates[config2].append(tr_rate)
    
    for config, rates in tr_rates.items():
        lifetime = 1/sum(rates)*1e9
        print(config.replace('.',' '), round(lifetime,2), 'ns')

if __name__ == "__main__":
    calc_lifetimes('tr_test.csv')