import pandas as pd
import requests

def get_radii():
    url = "https://www-nds.iaea.org/radii/"

    r = requests.get(url)
    df_list = pd.read_html(r.text)
    radii_df = df_list[0].drop([0, 1, 2])

    radii_df.columns = [ "Z", "Elem.", "Mass", "n", "R_av(fm)", "ΔR_av(fm)", "R_av,p(fm)", "ΔR_av,p(fm)"]
    radii_df.reset_index(drop=True, inplace=True)
    radii_df['Z'].fillna(method='ffill',inplace=True)
    radii_df['Elem.'].fillna(method='ffill',inplace=True)

    return radii_df

def get_extra_radii():
    rnuc_data = [
        [89,'Ac',227,5.7350],
        [91,'Pa',231,5.8000],
        [93,'Np',237,5.8500]
    ]
    rnuc_df = pd.DataFrame(rnuc_data, columns = ['Z','Elem.','Mass','R_av(fm)'])

    return rnuc_df

def get_periodic_table():
    url = "https://www.webelements.com/compounds.html"

    r = requests.get(url)
    df_list = pd.read_html(r.text)
    df = df_list[0]

    ptable = []

    # Append elements to a list ptable
    for index, row in df.iterrows():
        for col in df:
            if isinstance(row[col], str) and len(row[col]) > 12:
                element = row[col].split()
                # Remove radioactivity designation
                if len(element) == 5:
                    element.pop(2)
                ptable.append(element)

    # Create new dataframe for periodic table
    ptable_df = pd.DataFrame(ptable, columns = ['Z', 'Symbol', 'Mass', 'Name'])
    ptable_df.Z = ptable_df.Z.astype(int)
    ptable_df.Mass = ptable_df.Mass.astype(float)
    ptable_df.sort_values(by=['Z'], ascending=True, inplace=True)
    ptable_df.reset_index(drop=True, inplace=True)

    return ptable_df

def save_to_csv(df, filename):
    df.to_csv(filename)

if __name__ == "__main__":
    radii_df = get_radii()
    save_to_csv(radii_df, 'radii.csv')

    ptable_df = get_periodic_table()
    save_to_csv(ptable_df, 'ptable.csv')