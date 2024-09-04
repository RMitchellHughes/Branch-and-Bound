import pandas as pd
from BB_confound_tf import BB_confound_tf

def main():
    # Read in data and divide between x, y, and s
    dat = pd.read_csv('/Users/mitchell/Documents/NHANES Project/NHANES_07_12.csv')

    x = dat['SD.level']
    y = dat['bmi']
    s = dat.drop(['SEQN', 'SD.level', 'bmi'], axis = 1)

    # Run the BB algorithm
    bounds = BB_confound_tf(x, y, s)
    print(bounds)

if __name__ == '__main__':
    main()