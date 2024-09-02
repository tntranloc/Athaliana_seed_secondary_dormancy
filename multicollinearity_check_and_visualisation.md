```python
import os
import pandas as pd
import matplotlib as mpl
import seaborn as sns
```


```python
import patsy
import statsmodels
```


```python
from patsy import dmatrices
from statsmodels.stats.outliers_influence import variance_inflation_factor
```


```python
import numpy as np
```


```python
os.chdir('/working/dir')
```


```python
df = pd.read_csv("/your/dir/yourdata.csv")
```


```python

#find design matrix for regression model using 'rating' as response variable 
y, X = dmatrices('rate ~ bio3 + bio10+ bio9 + bio18', data=df, return_type='dataframe')

#create DataFrame to hold VIF values
vif_df = pd.DataFrame()
vif_df['variable'] = X.columns 

#calculate VIF for each predictor variable 
vif_df['VIF'] = [variance_inflation_factor(X.values, i) for i in range(X.shape[1])]

#view VIF for each predictor variable 
print(vif_df)



#In VIF (Variance Inflation Factor) analysis, 
#"const" typically refers to the intercept term in a regression model. 
#The intercept term represents the expected value of the dependent variable 
#when all the predictor variables are set to zero. 
#In VIF analysis, the intercept term is often included as a predictor variable, 
#which is why it appears as "const" in the results. 
#The VIF for "const" indicates the degree of multicollinearity 
#among the other predictor variables when controlling for the intercept term.
```

        variable         VIF
    0  Intercept  149.510707
    1       bio3    2.152819
    2      bio10    1.517482
    3       bio9    3.547181
    4      bio18    1.774200



```python
#visualising
#simple and fast option

# load the data
df1 = df[['bio9','bio10','bio18']]
# calculate the correlation matrix
corr = df1.corr(method="kendall")

#change lables if needed
#labels = {
#'bio9':'temperature of driest quarter',  
#'bio1':'annual temperature',
#'bio11':'temperature of coldest quarter',
#'bio6':'temperature of coldest month',
#'bio3':'isothermality'}

#corr = corr.rename(labels)

#make heatmap
myfig= sns.heatmap(corr, cmap="Blues", annot=True)

#figure = myfig.get_figure()    
#figure.savefig('multicorr_bio91812.png', dpi=800)
```


    
![png](output_7_0.png)
    

