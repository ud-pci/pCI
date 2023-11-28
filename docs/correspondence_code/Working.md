# Working

*Note : Definitions that performs the step explained is given in code blocks*

**The complete code works in four steps :**

## Correct Configurations

Before finding correspondance the configurations in UD data needs to be corrected. This is done in three steps

1. **Starting Configurations** - First configuration of each Term apearing in the NIST data is stored in a list as a starting point

```python
StartingConfigs # Starting Configurations
```

2. **Counting Configurations** - Every specific Term is counted to its multiplicity to make sure all conigurations are present

3. **Correcting Configurations** - In case a discrepancy is present, configuration will be corrected by subtracting the difference between discrepancy configuration and the previous configuration

|  Config | Term | J |     |  Config |
|    -    |   -  | - |  -  |    -    |
|  4d.5p  |  3D  | 1 | ==> |  4d.5p  |
|  4d.5p  |  3D  | 2 | ==> |  4d.5p  |
|**4d.7p**|  3D  | 3 | ==> |**4d.5p**|

```python
Correct_Config # Counting Configurations and Correcting Configurations
Corrected_Config # Calling "Correct_Config" and Storing Corrected Configuration
```

*Note : Some corrections are also need to be added manually which are not possibe using the above code. For the cases when the configuration is missidenified, such as if "4d.5p" is identified as "5s.9p"*

The corrected data is then used to find the Correspondance.

## Find Correspondance

Correspondance between a configuration in UD data is made with configurations present in the NIST data in five passes

Select a j<sup>th</sup> configuration in NIST data

1. **First Pass** - All configurations in UD with the J values equal to the J value of j<sup>th</sup> configuration in NIST data are selected in the first loop, ie US configurations with J<sub>UD</sub> = J<sub>NIST</sub>

2. **Second Pass** - From the configurations from the first pass those are selected with term same as term of j<sup>th</sup> configuration in NIST data, ie UD configurations with Term<sub>UD</sub> = Term<sub>NIST</sub>

3. **Third Pass** - In next step those configurations either "config1" or "config2" are selected which are same as of j<sup>th</sup> configuration in NIST data, ie Config<sub>UD</sub> = Config<sub>NIST</sub>

4. **Fourth Pass** - For the configurations for which term does't match in second pass, the correspondance is checked by adding and subtracting 1 from the integer value in term in UD data.

For example if term of UD configuration is "2D" and it failed at second pass then try checking with "1D" and "3D" and select this configuration if matched with NIST term of j<sup>th</sup> configuration

```python
FindJthAll # Above four passes are performed and returns j1,j2,dE1,dE2

# j1 - list of indexes of UD data for which config1 matched with jth NIST configuration
# j2 - list of indexes of UD data for which config2 matched with jth NIST configuration
# dE1 - List of percentage differencesin energies for configurations in j1 with the jth NIST configuration
# dE2 - List of percentage differencesin energies for configurations in j2 with the jth NIST configuration
```

5. **Fifth Pass** - Between configurations present in j1 and j2 the final selecten is made as follows
- Choose the configurations in j1 with minimum percentage difference if j1 is not empty and is the final configuration
  
- Select the configurations in j2 with minimum percentage difference if j2 is not empty and j1 is empty.
  - Check if Config1 corresponding the selected configuration is being used by some other NIST configuration with percentage difference less than the selected configuration.
  - If yes, leave the configuration for a better match
  - If no, choose this configuration as the final configuration

```python
ChooseJth # Choose between j1 and j2 to find the best correspondance
FindJth # Call "FindJthAll" and "ChooseJth" and returns index of final configuration surviving all five passes
```

## Combine Everything

1. Correct Configurations
2. Setup empty list (let us call it `"data_csv"`), to fill data
3. Setup a loop over all NIST states and call "FindJth" to get the corresponding configurations from corrected configurations of UD Data and add them to `"data_csv"`. Correspondance is shown by writing NIST configuration first and then corresponding UD configuration next to it.
4. NIST configurations with no corresponding UD configuration are then added to `"data_csv"` writing NIST configuration first and `"-"` next to mark it.
5. UD configurations with no corresponding NIST configuration are then added to `"data_csv"` writing `"-"` first to mark and then UD configuration next to `"-"`.

```python
MainCode # Combines all the steps 
```


## Ordering and Export

**Term Ordered** - The data is sorted by making the configurations ordered according to terms.

```python
FinalSortT # Sorting Data accoring to Terms
```

**Energy Ordered** - The data is sorted by making the configurations ordered according to NIST energies and using UD energies when NIST energy is not present for a configuration.

```python
FinalSortE # Sorting Data accoring to Energies
```


**Missing Configurations** - Once all steps are completed all the configurations and terms present in the data are then varified with configuration and terms calculated using online term calculator `http://umop.net/spectra/term_calc.php`. The configurations missing in both NIST and UD data are then added at the end of the data marked with "~" as "~3d2"

```python
Missing_Levels # Finding Missing Levels
```

**Export** - The final data is then exported in a txt format

```python
ConvertToTXT # Creating a txt file in DATA_Output
```