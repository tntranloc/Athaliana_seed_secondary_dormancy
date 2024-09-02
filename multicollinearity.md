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
os.chdir('/Users/nhutran/Documents/PhD/dormancy_redo')
```


```python
df = pd.read_csv("/Users/nhutran/Documents/PhD/dormancy_redo/assays/assay_sdorm02_working_bio.csv")
```


```python

#find design matrix for regression model using 'rating' as response variable 
y, X = dmatrices('rate ~ bio3 + bio10+ bio9 +bio18', data=df, return_type='dataframe')

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
    



```python
#make it fancy? 
#ver1
# labels need changing?
labels = {
'bio9':'temperature of driest quarter',  
'bio1':'annual temperature',
'bio11':'temperature of coldest quarter',
'bio6':'temperature of coldest month',
'bio3':'isothermality'}

corr = corr.rename(labels)

# remove the top right triange - duplicate information
mask = np.zeros_like(corr, dtype=np.bool)
mask[np.triu_indices_from(mask)] = True

# Colors
cmap = sns.diverging_palette(500, 10, as_cmap=True)

# uncomment this if you want only the lower triangle matrix 
# ans=sns.heatmap(corr, mask=mask,  linewidths=1, cmap=cmap, center=0)

ans=sns.heatmap(corr,  linewidths=1, cmap=cmap, center=0)

#save image 
#figure = ans.get_figure()    
#figure.savefig('correlations.png', dpi=800)
```

    /var/folders/_n/6f9lbyk959q4h617mnpxbxx40000gn/T/ipykernel_84123/1355055706.py:14: DeprecationWarning: `np.bool` is a deprecated alias for the builtin `bool`. To silence this warning, use `bool` by itself. Doing this will not modify any behavior and is safe. If you specifically wanted the numpy scalar type, use `np.bool_` here.
    Deprecated in NumPy 1.20; for more details and guidance: https://numpy.org/devdocs/release/1.20.0-notes.html#deprecations
      mask = np.zeros_like(corr, dtype=np.bool)



    
![png](output_8_1.png)
    



```python
#ver 2, a table 
cmap = cmap=sns.diverging_palette(5, 250, as_cmap=True)

def magnify():
    return [dict(selector="th",
                 props=[("font-size", "7pt")]),
            dict(selector="td",
                 props=[('padding', "0em 0em")]),
            dict(selector="th:hover",
                 props=[("font-size", "12pt")]),
            dict(selector="tr:hover td:hover",
                 props=[('max-width', '200px'),
                        ('font-size', '12pt')])
]

corr.style.background_gradient(cmap, axis=1)\
    .set_properties(**{'max-width': '80px', 'font-size': '10pt'})\
    .set_caption("Hover to magify")\
    .set_precision(2)\
    .set_table_styles(magnify())
```

    /var/folders/_n/6f9lbyk959q4h617mnpxbxx40000gn/T/ipykernel_84123/3442623433.py:16: FutureWarning: this method is deprecated in favour of `Styler.format(precision=..)`
      corr.style.background_gradient(cmap, axis=1)\





<style type="text/css">
#T_c06ea_ th {
  font-size: 7pt;
}
#T_c06ea_ td {
  padding: 0em 0em;
}
#T_c06ea_ th:hover {
  font-size: 12pt;
}
#T_c06ea_ tr:hover td:hover {
  max-width: 200px;
  font-size: 12pt;
}
#T_c06ea_row0_col0, #T_c06ea_row1_col1, #T_c06ea_row2_col2, #T_c06ea_row3_col3, #T_c06ea_row4_col4 {
  background-color: #4479bb;
  color: #f1f1f1;
  max-width: 80px;
  font-size: 10pt;
}
#T_c06ea_row0_col1, #T_c06ea_row1_col2, #T_c06ea_row2_col1, #T_c06ea_row3_col2, #T_c06ea_row4_col1 {
  background-color: #d73c5b;
  color: #f1f1f1;
  max-width: 80px;
  font-size: 10pt;
}
#T_c06ea_row0_col2 {
  background-color: #cfdae7;
  color: #000000;
  max-width: 80px;
  font-size: 10pt;
}
#T_c06ea_row0_col3 {
  background-color: #edcdd3;
  color: #000000;
  max-width: 80px;
  font-size: 10pt;
}
#T_c06ea_row0_col4 {
  background-color: #8aaad1;
  color: #f1f1f1;
  max-width: 80px;
  font-size: 10pt;
}
#T_c06ea_row1_col0 {
  background-color: #de6c83;
  color: #f1f1f1;
  max-width: 80px;
  font-size: 10pt;
}
#T_c06ea_row1_col3 {
  background-color: #f0dee2;
  color: #000000;
  max-width: 80px;
  font-size: 10pt;
}
#T_c06ea_row1_col4 {
  background-color: #e492a3;
  color: #f1f1f1;
  max-width: 80px;
  font-size: 10pt;
}
#T_c06ea_row2_col0 {
  background-color: #bccce1;
  color: #000000;
  max-width: 80px;
  font-size: 10pt;
}
#T_c06ea_row2_col3 {
  background-color: #efd7dc;
  color: #000000;
  max-width: 80px;
  font-size: 10pt;
}
#T_c06ea_row2_col4 {
  background-color: #5585c0;
  color: #f1f1f1;
  max-width: 80px;
  font-size: 10pt;
}
#T_c06ea_row3_col0 {
  background-color: #dc5c76;
  color: #f1f1f1;
  max-width: 80px;
  font-size: 10pt;
}
#T_c06ea_row3_col1 {
  background-color: #d94764;
  color: #f1f1f1;
  max-width: 80px;
  font-size: 10pt;
}
#T_c06ea_row3_col4 {
  background-color: #e7a8b5;
  color: #000000;
  max-width: 80px;
  font-size: 10pt;
}
#T_c06ea_row4_col0 {
  background-color: #95b1d5;
  color: #000000;
  max-width: 80px;
  font-size: 10pt;
}
#T_c06ea_row4_col2 {
  background-color: #5c8ac2;
  color: #f1f1f1;
  max-width: 80px;
  font-size: 10pt;
}
#T_c06ea_row4_col3 {
  background-color: #f1e7e9;
  color: #000000;
  max-width: 80px;
  font-size: 10pt;
}
</style>
<table id="T_c06ea_">
  <caption>Hover to magify</caption>
  <thead>
    <tr>
      <th class="blank level0" >&nbsp;</th>
      <th class="col_heading level0 col0" >bio1</th>
      <th class="col_heading level0 col1" >bio3</th>
      <th class="col_heading level0 col2" >bio6</th>
      <th class="col_heading level0 col3" >bio9</th>
      <th class="col_heading level0 col4" >bio11</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th id="T_c06ea_level0_row0" class="row_heading level0 row0" >bio1</th>
      <td id="T_c06ea_row0_col0" class="data row0 col0" >1.00</td>
      <td id="T_c06ea_row0_col1" class="data row0 col1" >0.55</td>
      <td id="T_c06ea_row0_col2" class="data row0 col2" >0.82</td>
      <td id="T_c06ea_row0_col3" class="data row0 col3" >0.73</td>
      <td id="T_c06ea_row0_col4" class="data row0 col4" >0.91</td>
    </tr>
    <tr>
      <th id="T_c06ea_level0_row1" class="row_heading level0 row1" >bio3</th>
      <td id="T_c06ea_row1_col0" class="data row1 col0" >0.55</td>
      <td id="T_c06ea_row1_col1" class="data row1 col1" >1.00</td>
      <td id="T_c06ea_row1_col2" class="data row1 col2" >0.48</td>
      <td id="T_c06ea_row1_col3" class="data row1 col3" >0.71</td>
      <td id="T_c06ea_row1_col4" class="data row1 col4" >0.61</td>
    </tr>
    <tr>
      <th id="T_c06ea_level0_row2" class="row_heading level0 row2" >bio6</th>
      <td id="T_c06ea_row2_col0" class="data row2 col0" >0.82</td>
      <td id="T_c06ea_row2_col1" class="data row2 col1" >0.48</td>
      <td id="T_c06ea_row2_col2" class="data row2 col2" >1.00</td>
      <td id="T_c06ea_row2_col3" class="data row2 col3" >0.70</td>
      <td id="T_c06ea_row2_col4" class="data row2 col4" >0.97</td>
    </tr>
    <tr>
      <th id="T_c06ea_level0_row3" class="row_heading level0 row3" >bio9</th>
      <td id="T_c06ea_row3_col0" class="data row3 col0" >0.73</td>
      <td id="T_c06ea_row3_col1" class="data row3 col1" >0.71</td>
      <td id="T_c06ea_row3_col2" class="data row3 col2" >0.70</td>
      <td id="T_c06ea_row3_col3" class="data row3 col3" >1.00</td>
      <td id="T_c06ea_row3_col4" class="data row3 col4" >0.79</td>
    </tr>
    <tr>
      <th id="T_c06ea_level0_row4" class="row_heading level0 row4" >bio11</th>
      <td id="T_c06ea_row4_col0" class="data row4 col0" >0.91</td>
      <td id="T_c06ea_row4_col1" class="data row4 col1" >0.61</td>
      <td id="T_c06ea_row4_col2" class="data row4 col2" >0.97</td>
      <td id="T_c06ea_row4_col3" class="data row4 col3" >0.79</td>
      <td id="T_c06ea_row4_col4" class="data row4 col4" >1.00</td>
    </tr>
  </tbody>
</table>



