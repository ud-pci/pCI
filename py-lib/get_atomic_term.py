import re
import requests
import atomic_term_symbol

def strip_text(str, txt_to_remove):
    '''
    This function strips a string of unwanted text
    '''
    for txt in txt_to_remove:
        str = str.replace(txt, '')
        
    return str

def scrape_term(configuration):
    '''
    This function calculates all possible terms for a given electronic configuration
    '''
    all_terms = atomic_term_symbol.calc_term_symbols(configuration)
    terms = []
    for term in all_terms:
        if ',' in term:
            nl = term[:2]
            J_list = term[2:].split(',')
            
            for J in J_list:
                terms.append(nl + J)
        else:
            terms.append(term)
    
    return terms
    

if __name__ == '__main__':
    configuration = input('Input configuration: ')
    terms = scrape_term(configuration)
    print(len(terms), 'terms for ' + configuration + ':', terms)