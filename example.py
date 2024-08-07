# An example application of our algorithm (see section 4)

import pandas as pd
from BB_confound import BB_confound

def main():
    # Read in data and divide between x, y, and s
    dat = pd.read_csv('/Users/mitchell/Documents/NHANES Project/NHANES_07_12.csv')

    x = dat['SD.level']
    y = dat['bmi']
    s = dat.drop(['SEQN', 'SD.level', 'bmi'], axis = 1)

    # Create interaction terms
    s_interact = pd.DataFrame()
    for var1 in range(len(s.columns)):
        for var2 in range(var1,len(s.columns)):
            if var1 == var2:
                continue
            var_name = f"{s.columns[var1]}_{s.columns[var2]}"
            s_interact[var_name] = s.iloc[:,var1] * s.iloc[:,var2]

    # Eliminate interaction between racial indicator variables
    s_interact = s_interact.drop(['white_black'], axis = 1)
    s = pd.concat([s, s_interact], axis = 1)

    # Run the BB algorithm
    bounds = BB_confound(x, y, s)
    print(bounds)

if __name__ == '__main__':
    main()

# Note: Because the data set is large, it may take a while for the algorithm to finish.