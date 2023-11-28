import re
import requests

def strip_text(str, txt_to_remove):
    '''
    This function strips a string of unwanted text
    '''
    for txt in txt_to_remove:
        str = str.replace(txt, '')
        
    return str

def scrape_term(configuration):
    '''
    This function scrapes all possible terms for a given electron configuration from online term calculator
    at http://umop.net/spectra/term_calc.php
    '''
    url = 'http://umop.net/spectra/term_calc.php?config=' + configuration.replace(' ', '.')
    r = requests.get(url)
    
    reg_str = "<strong>(.*?)</strong>"
    txt_to_remove = ['<sub>','</sub>','<sup>','</sup>','&nbsp;']
    scraped_terms = strip_text(re.findall(reg_str, r.text)[0].replace('&deg;',''), txt_to_remove).split()
    
    terms = []
    for term in scraped_terms:
        if ',' in term:
            nl = term[:2]
            J_list = term[2:].split(',')
            
            for J in J_list:
                terms.append(nl + J)
        else:
            terms.append(term)
    
    return terms
    

# if __name__ == '__main__':
#     configuration = input('Input configuration: ')
#     terms = scrape_term(configuration)
#     print(len(terms), 'terms for ' + configuration + ':', terms)