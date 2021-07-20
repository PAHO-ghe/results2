[Home](https://paho-ghe.github.io/PAHO/)
# COVID19 trend anormally detection 

## 1.Check COVID19 trend anormally in PAHO countries

This function (forked from github trendbreaker) implements an algorithm for epidemic time series analysis in aim to detect recent deviation from the trend followed by the data. Data is first partitioned into 'recent' data, using the last k observations as supplementary individuals, and older data used to fit the trend. Trend-fitting is done by fitting a series of user-specified models for the time series, with different methods for selecting best fit (see details, and the argument method). The prediction interval is then calculated for the best model, and every data point (including the training set and supplementary individuals) falling outside are classified as 'outliers'. The value of k can be fixed by the user, or automatically selected to minimise outliers in the training period and maximise and the detection of outliers in the recent period.

``` r
setwd("./")
GHE_PAHO<-read.csv("./GHE_PAHO.csv")
GHE_PAHO_NDeaths_AgeGroup     <- GHE_PAHO[GHE_PAHO$ghecause==0,-c(5,6,8,9)] 
GHE_PAHO_NDeaths              <- data.frame(xtabs(dths~iso3+year, GHE_PAHO_NDeaths_AgeGroup))
names(GHE_PAHO_NDeaths)[3]    <- "dths"

#save(GHE_PAHO_NDeaths, file = "./GHE_PAHO_NDeaths.RData")
#load(GHE_PAHO_NDeaths, file = "./GHE_PAHO_NDeaths.RData")
```

## 2.Prep the module

``` r
GHE_PAHO_NDeaths$year <- as.Date(GHE_PAHO_NDeaths$year, format = "%Y")
GHE_PAHO_NDeaths$dths <- as.integer(GHE_PAHO_NDeaths$dths)

first_date <- max(GHE_PAHO_NDeaths$year, na.rm = TRUE) - 10
pathways_recent <- GHE_PAHO_NDeaths %>%filter(year >= first_date)

# define candidate models
models <- list(
  regression = lm_model(dhts ~ year), 
  poisson_constant = glm_model(dhts ~ year, family = "poisson"),
  negbin_time = glm_nb_model(dhts ~ year) #negative binomial generalied linear models
)

# analyses on all data
counts_overall <- pathways_recent %>%
  group_by(year) %>%
  summarise(count = sum(dths))
# results with automated detection of 'k'
#res <- asmodee(counts_overall, models, max_k = 2, date_index = "year", method = evaluate_aic)
#res
```


> Impression : Dimension of data not applicable for this analysis's max\_k window size.

## 3. Run it on PAHO region COVID19 daily number of deaths

``` r
# get data
data("coronavirus")

Country_codes     <- readxl::read_excel("./Country_codes_ISO_3166.xlsx")
Country_codes<-Country_codes%>%select(nom.eng,Continent)%>%filter(Continent=="Americas")
df_paho<-merge(coronavirus, Country_codes, by.x="country", by.y="nom.eng")

df_paho<-df_paho%>%select(date, country,cases,type)%>%filter(type=="death")
#df_arg<-df_paho%>%select(date, country,cases,type)%>%filter(country=="Argentina") # Example with one of the large population country
#dim(df_arg)
df_arg<-df_paho
```

``` r
#prep data
df_arg<-df_paho
df_arg<-df_arg%>%arrange(date)
df_arg$cases <- as.integer(df_arg$cases)
df_arg$day<- difftime(df_arg$date,min(df_arg$date),  units = c("days"))
df_arg$day <- as.integer(df_arg$day)
df_arg$weekday <- weekdays(as.Date(df_arg$date))
df_arg$weekday<-recode(df_arg$weekday, Saturday="weekend", Sunday='weekend',Tuesday='rest_of_week',Wednesday='rest_of_week',Thursday='rest_of_week',Friday='rest_of_week',Monday='monday')
df_arg$weekday <- as.factor(df_arg$weekday)
df_arg<-as_tibble(df_arg%>%select(date,cases,day,weekday))%>%arrange(date)
names(df_arg)[2]    <- "count"

first_date <- max(df_arg$date, na.rm = TRUE) - 28
pathways_recent <- df_arg %>%filter(date >= first_date)

models <- list(
    regression = lm_model(count ~ day),
    poisson_constant = glm_model(count ~ day, family = "poisson"),
    negbin_time = glm_nb_model(count ~ day),
    negbin_time_weekday = glm_nb_model(count ~ day + weekday)
              )
counts_overall <- pathways_recent %>%group_by(date, day, weekday) %>%summarise(count = sum(count))
```

    ## `summarise()` has grouped output by 'date', 'day'. You can override using the `.groups` argument.

``` r
res_overall <- asmodee(counts_overall,models,"date",method = evaluate_aic)
res_overall$fitted_model
```

    ## NULL

``` r
res_overall$results
```
|	date        | day          | training	| estimate  | lower_ci | upper_ci	|	lower_pi  | upper_pi |	outlier	|	classification |
|	----------- |	-------------- | ------ | ---------- | -------- | -------- | --------- | --------- | --------- |-----------------|
|	2021-01-01 	|	 weekday count	|	TRUE	|	  2579.499 | 2189.247 | 3039.317	|	1642	|	3755	|	  FALSE 	|	        normal	|
|	2021-01-02 	|	rest_of_week  	|	TRUE	|	  1783.453 | 1504.439 | 2114.213	|	1120	|	2529	|	  FALSE 	|	        normal	|
|	2021-01-03 	|	     weekend  	|	TRUE	|	  1819.393 | 1542.320 | 2146.242	|	1143	|	2629	|	  FALSE 	|	        normal	|
|	2021-01-04 	|	     weekend  	|	TRUE	|	 1734.652 |1394.785 | 2157.336	  |	1106	|	2610	|	  FALSE 	|	        normal	|
|	2021-01-05 	|	      monday  	|	TRUE	|	 2793.798 | 2436.622 | 3203.332	|	1814	|	3969	|	  FALSE 	|	        normal	|
|	2021-01-06 	|	rest_of_week  	|	TRUE	|	 2850.099 |2501.180 | 3247.693	  |	1871	|	4081	|	  FALSE 	|	        normal	|
|	2021-01-07 	|	rest_of_week  	|	TRUE	|	  2907.535 | 2566.625 | 3293.726	  |	1900	|	4226	|	  FALSE 	|	        normal	|
|	2021-01-08 	|	rest_of_week  	|	TRUE	|	 2966.128 |2632.812 | 3341.641	  |	1942	|	4352	|	  FALSE 	|	        normal	|
|	2021-01-09 	|	rest_of_week  	|	TRUE	|	 2050.766 |1775.975 | 2368.075	  |	1307	|	2968	|	TRUE    	|	      increase	|
|	2021-01-10 	|	     weekend  	|	TRUE	|	 2092.093 |1815.262 | 2411.142	  |	1349	|	2988	|	  FALSE 	|	        normal	|
|	2021-01-11 	|	     weekend  	|	TRUE	|	 1994.651 |1633.681 | 2435.379	  |	1250	|	3036	|	  FALSE 	|	        normal	|
|	2021-01-12 	|	      monday  	|	TRUE	|	 3212.547 | 2901.043 | 3557.500	|	2062	|	4463	|	  FALSE 	|	        normal	|
|	2021-01-13 	|	rest_of_week  	|	TRUE	|	 3277.287 | 2967.750 | 3619.108  	|	2147	|	4750	|	  FALSE 	|	        normal	|
|	2021-01-14 	|	rest_of_week  	|	TRUE	|	  3343.331 | 3033.776 | 3684.471  	|	2220	|	4694	|	  FALSE 	|	        normal	|
|	2021-01-15 	|	rest_of_week  	|	TRUE	|	 3410.706 | 3098.863 | 3753.930  	|	2235	|	4814	|	  FALSE 	|	        normal	|
|	2021-01-16 	|	rest_of_week  	|	TRUE	|	 2358.145 | 2046.212 | 2717.631	  |	1505	|	3453	|	  FALSE 	|	        normal	|
|	2021-01-17 	|	     weekend  	|	TRUE	|	  2405.667 | 2083.459 | 2777.704	  |	1567	|	3480	|	  FALSE 	|	        normal	|
|	2021-01-18 	|	     weekend  	|	TRUE	|	 2293.619 | 1878.626 | 2800.285	  |	1452	|	3325	|	  FALSE 	|	        normal	|
|	2021-01-19 	|	      monday  	|	TRUE	|	  3694.060 | 3346.279 | 4077.987	  |	2406	|	5371	|	  FALSE 	|	        normal	|
|	2021-01-20 	|	rest_of_week  	|	TRUE	|	 3768.503 | 3404.568 | 4171.342	  |	2455	|	5461	|	  FALSE 	|	        normal	|
|	2021-01-21 	|	rest_of_week  	|	TRUE	|	 3844.447 | 3461.515 | 4269.740	  |	2552	|	5537	|	  FALSE 	|	        normal	|
|	2021-01-22 	|	rest_of_week  	|	TRUE	|	 3921.920 | 3517.241 | 4373.160	  |	2527	|	5540	|	  FALSE 	|	        normal	|
|	2021-01-23 	|	rest_of_week  	|	TRUE	|	 2711.596 | 2298.981 | 3198.266	  |	1782	|	3999	|	  FALSE 	|	        normal	|
|	2021-01-24 	|	     weekend  	|	TRUE	|	 2766.240 | 2333.832 | 3278.765	  |	1771	|	4039	|	  FALSE 	|	        normal	|
|	2021-01-25 	|	     weekend  	|	TRUE	|	 2637.399 | 2120.909 | 3279.666	  |	1607	|	3950	|	  FALSE 	|	        normal	|
|	2021-01-26 	|	      monday  	|	TRUE	|	 4247.745 | 3730.885 | 4836.208	  |	2794	|	6229	|	  FALSE 	|	        normal	|
|	2021-01-27 	|	rest_of_week  	|	TRUE	|	 4333.346 | 3782.713 | 4964.132	  |	2870	|	6294	|	  FALSE 	|	        normal	|
|	2021-01-28 	|	rest_of_week  	|	TRUE	|	 4420.672 | 3834.169 | 5096.891	  |	2809	|	6532	|	  FALSE 	|	        normal	|
|	2021-01-29 	|	rest_of_week  	|	TRUE	|	  4509.758 | 3885.361 | 5234.497	  |	2940	|	6532	|	  FALSE 	|	        normal	|



![unnamed-chunk-4-1](https://user-images.githubusercontent.com/81782228/126196259-a16feb74-7c61-42ce-baf4-543235493a1d.png)

## 4. PAHO region COVID19 death trend in last 100 days since the last reported date (2021-01-29)

``` r
# try with wider range
df_arg<-df_paho
df_arg<-df_arg%>%arrange(date)
df_arg$cases <- as.integer(df_arg$cases)
df_arg$day<- difftime(df_arg$date,min(df_arg$date),  units = c("days"))
df_arg$day <- as.integer(df_arg$day)
df_arg$weekday <- weekdays(as.Date(df_arg$date))
df_arg$weekday<-recode(df_arg$weekday, Saturday="weekend", Sunday='weekend',Tuesday='rest_of_week',Wednesday='rest_of_week',Thursday='rest_of_week',Friday='rest_of_week',Monday='monday')
df_arg$weekday <- as.factor(df_arg$weekday)
df_arg<-as_tibble(df_arg%>%select(date,cases,day,weekday))%>%arrange(date)
names(df_arg)[2]    <- "count"

first_date <- max(df_arg$date, na.rm = TRUE) - 100
pathways_recent <- df_arg %>%filter(date >= first_date)

models <- list(
  regression = lm_model(count ~ day),
  poisson_constant = glm_model(count ~ day, family = "poisson"),
  negbin_time = glm_nb_model(count ~ day),
  negbin_time_weekday = glm_nb_model(count ~ day + weekday)
)
counts_overall <- pathways_recent %>%group_by(date, day, weekday) %>%summarise(count = sum(count))
```

    ## `summarise()` has grouped output by 'date', 'day'. You can override using the `.groups` argument.

``` r
res_overall <- asmodee(counts_overall,models,"date",method = evaluate_aic)
res_overall$fitted_model
```

    ## NULL

``` r
res_overall$results
```

|	      date 	|	     weekday 	|	count |	training 	|	estimate 	|	 lower_ci |	upper_ci	|	lower_pi |	upper_pi |	outlier	|	classification	|
| ----------- | ------------ | ------ | ---------- | ---------- | ---------- | ---------- | -------- | ---------- | -------- | --------------- |
|	2020-10-21 	|	rest_of_week 	|	1953	|	    TRUE 	|	1599.798	|	1466.4433	|	1745.279	|	1023	  |	2351	|	FALSE 	|	normal	|
|	2020-10-22 	|	rest_of_week 	|	1452	|	    TRUE 	|	1613.056	|	1480.2166	|	1757.818	|	1024	  |	2343	|	FALSE 	|	normal	|
|	2020-10-23 	|	rest_of_week 	|	2183	|	    TRUE 	|	1626.425	|	1494.1069	|	1770.461	|	1033	  |	2346	|	FALSE 	|	normal	|
|	2020-10-24 	|	     weekend 	|	1554	|	    TRUE 	|	1095.216	|	991.6890	|	1209.551	|	705	    |	1564	|	FALSE 	|	normal	|
|	2020-10-25 	|	     weekend 	|	1137	|	    TRUE 	|	1104.293	|	1000.7829	|	1218.509	|	703	    |	1607	|	FALSE 	|	normal	|
|	2020-10-26 	|	      monday 	|	1316	|	    TRUE 	|	1069.721	|	944.0906	|	1212.068	|	687	    |	1544	|	FALSE 	|	normal	|
|	2020-10-27 	|	rest_of_week 	|	2028	|	    TRUE 	|	1681.016	|	1550.8461	|	1822.112	|	1086	  |	2427	|	FALSE 	|	normal	|
|	2020-10-28 	|	rest_of_week 	|	1744	|	    TRUE 	|	1694.948	|	1565.3268	|	1835.302	|	1093	  |	2437	|	FALSE 	|	normal	|
|	2020-10-29 	|	rest_of_week 	|	1805	|	    TRUE 	|	1708.995	|	1579.9264	|	1848.607	|	1083	  |	2536	|	FALSE 	|	normal	|
|	2020-10-30 	|	rest_of_week 	|	1802	|	    TRUE 	|	1723.159	|	1594.6450	|	1862.029	|	1087	  |	2489	|	FALSE 	|	normal	|
|	2020-10-31 	|	     weekend 	|	1511	|	    TRUE 	|	1160.355	|	1056.8557	|	1273.991	|	736	    |	1666	|	FALSE 	|	normal	|
|	2020-11-01 	|	     weekend 	|	876	  |	    TRUE 	|	1169.972	|	1066.4534	|	1283.539	|	758	    |	1716	|	FALSE 	|	normal	|
|	2020-11-02 	|	      monday 	|	1170	|	    TRUE 	|	1133.344	|	1004.8104	|	1278.319	|	713	    |	1682	|	FALSE 	|	normal	|
|	2020-11-03 	|	rest_of_week 	|	1523	|	    TRUE 	|	1780.996	|	1654.7102	|	1916.921	|	1126	  |	2586	|	FALSE 	|	normal	|
|	2020-11-04 	|	rest_of_week 	|	2226	|	    TRUE 	|	1795.757	|	1670.0240	|	1930.956	|	1170	  |	2671	|	FALSE 	|	normal	|
|	2020-11-05 	|	rest_of_week 	|	1242	|	    TRUE 	|	1810.639	|	1685.4565	|	1945.120	|	1160	  |	2613	|	FALSE 	|	normal	|
|	2020-11-06 	|	rest_of_week 	|	2243	|	    TRUE 	|	1825.645	|	1701.0076	|	1959.416	|	1181	  |	2604	|	FALSE 	|	normal	|
|	2020-11-07 	|	     weekend 	|	1405	|	    TRUE 	|	1229.369	|	1125.5555	|	1342.757	|	792	    |	1743	|	FALSE 	|	normal	|
|	2020-11-08 	|	     weekend 	|	950	  |	    TRUE 	|	1239.557	|	1135.6579	|	1352.963	|	795	    |	1803	|	FALSE 	|	normal	|
|	2020-11-09 	|	      monday 	|	1213	|	    TRUE 	|	1200.751	|	1068.7712	|	1349.028	|	761	    |	1761	|	FALSE 	|	normal	|
|	2020-11-10 	|	rest_of_week 	|	1458	|	    TRUE 	|	1886.923	|	1764.3886	|	2017.968	|	1213  	|	2699	|	FALSE 	|	normal	|
|	2020-11-11 	|	rest_of_week 	|	1968	|	    TRUE 	|	1902.562	|	1780.5256	|	2032.962	|	1215		|	2712	|	FALSE 	|	normal	|
|	2020-11-12 	|	rest_of_week 	|	2235	|	    TRUE 	|	1918.329	|	1796.7780	|	2048.104	|	1238		|	2744	|	FALSE 	|	normal	|
|	2020-11-13 	|	rest_of_week 	|	1728	|	    TRUE 	|	1934.228	|	1813.1450	|	2063.397	|	1246		|	2727	|	FALSE 	|	normal	|
|	2020-11-14 	|	     weekend 	|	2228	|	    TRUE 	|	1302.487	|	1197.7750	|	1416.353	|	832			|	1871	|	 TRUE 	|	increase	|
|	2020-11-15 	|	     weekend 	|	995	  |	    TRUE 	|	1313.282	|	1208.3765	|	1427.294	|	829			|	1884	|	FALSE 	|	normal	|
|	2020-11-16 	|	      monday 	|	1281	|	    TRUE 	|	1272.167	|	1136.0234	|	1424.626	|	801			|	1845	|	FALSE 	|	normal	|
|	2020-11-17 	|	rest_of_week 	|	1584	|	    TRUE 	|	1999.150	|	1879.7404	|	2126.146	|	1295		|	2905	|	FALSE 	|	normal	|
|	2020-11-18 	|	rest_of_week 	|	1983	|	    TRUE 	|	2015.719	|	1896.6657	|	2142.245	|	1254		|	2911	|	FALSE 	|	normal	|
|	2020-11-19 	|	rest_of_week 	|	1813	|	    TRUE 	|	2032.424	|	1913.6991	|	2158.515	|	1282		|	2901	|	FALSE 	|	normal	|
|	2020-11-20 	|	rest_of_week 	|	2037	|	    TRUE 	|	2049.268	|	1930.8392	|	2174.962	|	1341		|	2942	|	FALSE 	|	normal	|
|	2020-11-21 	|	     weekend 	|	1501	|	    TRUE 	|	1379.954	|	1273.4566	|	1495.358	|	883			|	1982	|	FALSE 	|	normal	|
|	2020-11-22 	|	     weekend 	|	969	  |	    TRUE 	|	1391.391	|	1284.5457	|	1507.123	|	887			|	1995	|	FALSE 	|	normal	|
|	2020-11-23 	|	      monday 	|	1150	|	    TRUE 	|	1347.830	|	1206.6061	|	1505.584	|	849			|	1939	|	FALSE 	|	normal	|
|	2020-11-24 	|	rest_of_week 	|	2235	|	    TRUE 	|	2118.052	|	2000.4376	|	2242.582	|	1354		|	3121	|	FALSE 	|	normal	|
|	2020-11-25 	|	rest_of_week 	|	2191	|	    TRUE 	|	2135.606	|	2018.0885	|	2259.967	|	1356		|	3094	|	FALSE 	|	normal	|
|	2020-11-26 	|	rest_of_week 	|	2013	|	    TRUE 	|	2153.305	|	2035.8366	|	2277.552	|	1375		|	3112	|	FALSE 	|	normal	|
|	2020-11-27 	|	rest_of_week 	|	1965	|	    TRUE 	|	2171.151	|	2053.6801	|	2295.342	|	1391		|	3110	|	FALSE 	|	normal	|
|	2020-11-28 	|	     weekend 	|	1745	|	    TRUE 	|	1462.028	|	1352.5078	|	1580.418	|	953			|	2129	|	FALSE 	|	normal	|
|	2020-11-29 	|	     weekend 	|	1067	|	    TRUE 	|	1474.145	|	1364.0697	|	1593.104	|	943			|	2136	|	FALSE 	|	normal	|
|	2020-11-30 	|	      monday 	|	1303	|	    TRUE 	|	1427.994	|	1280.5517	|	1592.413	|	932			|	2055	|	FALSE 	|	normal	|
|	2020-12-01 	|	rest_of_week 	|	2128	|	    TRUE 	|	2244.026	|	2125.9716	|	2368.636	|	1396		|	3197	|	FALSE 	|	normal	|
|	2020-12-02 	|	rest_of_week 	|	2158	|	    TRUE 	|	2262.624	|	2144.2648	|	2387.516	|	1476		|	3276	|	FALSE 	|	normal	|
|	2020-12-03 	|	rest_of_week 	|	2138	|	    TRUE 	|	2281.376	|	2162.6427	|	2406.627	|	1440		|	3328	|	FALSE 	|	normal	|
|	2020-12-04 	|	rest_of_week 	|	2063	|	    TRUE 	|	2300.283	|	2181.1036	|	2425.975	|	1474		|	3232	|	FALSE 	|	normal	|
|	2020-12-05 	|	     weekend 	|	1905	|	    TRUE 	|	1548.984	|	1434.8238	|	1672.228	|	984			|	2274	|	FALSE 	|	normal	|
|	2020-12-06 	|	     weekend 	|	1136	|	    TRUE 	|	1561.822	|	1446.8443	|	1685.937	|	996			|	2268	|	FALSE 	|	normal	|
|	2020-12-07 	|	      monday 	|	1364	|	    TRUE 	|	1512.926	|	1357.8933	|	1685.659	|	959			|	2204	|	FALSE 	|	normal	|
|	2020-12-08 	|	rest_of_week 	|	2147	|	    TRUE 	|	2377.492	|	2255.7458	|	2505.809	|	1549		|	3434	|	FALSE 	|	normal	|
|	2020-12-09 	|	rest_of_week 	|	2359	|	    TRUE 	|	2397.196	|	2274.5993	|	2526.401	|	1545		|	3471	|	FALSE 	|	normal	|
|	2020-12-10 	|	rest_of_week 	|	2278	|	    TRUE 	|	2417.063	|	2293.5275	|	2547.253	|	1530		|	3472	|	FALSE 	|	normal	|
|	2020-12-11 	|	rest_of_week 	|	1962	|	    TRUE 	|	2437.095	|	2312.5297	|	2568.370	|	1541		|	3526	|	FALSE 	|	normal	|
|	2020-12-12 	|	     weekend 	|	1935	|	    TRUE 	|	1641.112	|	1520.3201	|	1771.501	|	1075		|	2394	|	FALSE 	|	normal	|
|	2020-12-13 	|	     weekend 	|	1014	|	    TRUE 	|	1654.713	|	1532.7908	|	1786.333	|	1045		|	2380	|	 TRUE 	|	decrease	|
|	2020-12-14 	|	      monday 	|	1650	|	    TRUE 	|	1602.909	|	1438.6739	|	1785.893	|	1028		|	2316	|	FALSE 	|	normal	|
|	2020-12-15 	|	rest_of_week 	|	2464	|	    TRUE 	|	2518.896	|	2389.2653	|	2655.561	|	1645		|	3592	|	FALSE 	|	normal	|
|	2020-12-16 	|	rest_of_week 	|	2246	|	    TRUE 	|	2539.772	|	2408.6291	|	2678.056	|	1639		|	3659	|	FALSE 	|	normal	|
|	2020-12-17 	|	rest_of_week 	|	2672	|	    TRUE 	|	2560.821	|	2428.0648	|	2700.836	|	1637		|	3650	|	FALSE 	|	normal	|
|	2020-12-18 	|	rest_of_week 	|	2236	|	    TRUE 	|	2582.044	|	2447.5728	|	2723.904	|	1629		|	3676	|	FALSE 	|	normal	|
|	2020-12-19 	|	     weekend 	|	1935	|	    TRUE 	|	1738.719	|	1608.9668	|	1878.935	|	1114		|	2524	|	FALSE 	|	normal	|
|	2020-12-20 	|	     weekend 	|	1417	|	    TRUE 	|	1753.129	|	1621.8896	|	1894.988	|	1091		|	2497	|	FALSE 	|	normal	|
|	2020-12-21 	|	      monday 	|	1679	|	    TRUE 	|	1698.244	|	1522.9552	|	1893.708	|	1092		|	2458	|	FALSE 	|	normal	|
|	2020-12-22 	|	rest_of_week 	|	2731	|	    TRUE 	|	2668.711	|	2526.3368	|	2819.108	|	1688		|	3836	|	FALSE 	|	normal	|
|	2020-12-23 	|	rest_of_week 	|	2413	|	    TRUE 	|	2690.828	|	2546.2149	|	2843.655	|	1752		|	3822	|	FALSE 	|	normal	|
|	2020-12-24 	|	rest_of_week 	|	2264	|	    TRUE 	|	2713.129	|	2566.1700	|	2868.504	|	1781		|	3949	|	FALSE 	|	normal	|
|	2020-12-25 	|	rest_of_week 	|	1814	|	    TRUE 	|	2735.614	|	2586.2033	|	2893.657	|	1756		|	3919	|	FALSE 	|	normal	|
|	2020-12-26 	|	     weekend 	|	1102	|	    TRUE 	|	1842.132	|	1700.8125	|	1995.193	|	1158		|	2679	|	 TRUE 	|	decrease	|
|	2020-12-27 	|	     weekend 	|	1371	|	    TRUE 	|	1857.399	|	1714.2013	|	2012.558	|	1182		|	2676	|	FALSE 	|	normal	|
|	2020-12-28 	|	      monday 	|	1695	|	    TRUE 	|	1799.249	|	1610.8255	|	2009.713	|	1164		|	2573	|	FALSE 	|	normal	|
|	2020-12-29 	|	rest_of_week 	|	2939	|	    TRUE 	|	2827.436	|	2667.1481	|	2997.356	|	1830		|	4002	|	FALSE 	|	normal	|
|	2020-12-30 	|	rest_of_week 	|	3018	|	    TRUE 	|	2850.868	|	2687.5952	|	3024.061	|	1790		|	4054	|	FALSE 	|	normal	|
|	2020-12-31 	|	rest_of_week 	|	2839	|	    TRUE 	|	2874.496	|	2708.1302	|	3051.081	|	1864		|	4061	|	FALSE 	|	normal	|
|	2021-01-01 	|	rest_of_week 	|	1727	|	    TRUE 	|	2898.318	|	2728.7549	|	3078.419	|	1822		|	4173	|	 TRUE 	|	decrease	|
|	2021-01-02 	|	     weekend 	|	1249	|	    TRUE 	|	1951.695	|	1795.9908	|	2120.897	|	1290		|	2789	|	 TRUE 	|	decrease	|
|	2021-01-03 	|	     weekend 	|	1297	|	    TRUE 	|	1967.870	|	1809.8711	|	2139.661	|	1237	|	2815	|	FALSE 	|	normal	|
|	2021-01-04 	|	      monday 	|	1910	|	    TRUE 	|	1906.261	|	1702.4039	|	2134.530	|	1179	|	2773	|	FALSE 	|	normal	|
|	2021-01-05 	|	rest_of_week 	|	2949	|	    TRUE 	|	2995.601	|	2812.1886	|	3190.975	|	1945	|	4352	|	FALSE |	normal	|
|	2021-01-06 	|	rest_of_week 	|	3353	|	    TRUE 	|	3020.427	|	2833.2905	|	3219.924	|	1921	|	4303	|	FALSE 	|	normal	|
|	2021-01-07 	|	rest_of_week 	|	3475	|	    TRUE 	|	3045.460	|	2854.4939	|	3249.201	|	1942	|	4354	|	FALSE	|	normal	|
|	2021-01-08 	|	rest_of_week 	|	3042	|	    TRUE 	|	3070.699	|	2875.8007	|	3278.807	|	1930	|	4462	|	FALSE |	normal	|
|	2021-01-09 	|	     weekend 	|	3146	|	    TRUE 	|	2067.774	|	1894.7098	|	2256.646	|	1321	|	2939	|	 TRUE |	increase	|
|	2021-01-10 	|	     weekend 	|	1662	|	    TRUE 	|	2084.911	|	1909.1164	|	2276.893	|	1334	|	3049	|	FALSE |	normal	|
|	2021-01-11 	|	      monday 	|	2213	|	    TRUE 	|	2019.639	|	1797.8427	|	2268.797	|	1263	|	2954	|	FALSE |	normal	|
|	2021-01-12 	|	rest_of_week 	|	3375	|	    TRUE 	|	3173.768	|	2962.1019	|	3400.559	|	2043	|	4662	|	FALSE |	normal	|
|	2021-01-13 	|	rest_of_week 	|	3438	|	    TRUE 	|	3200.071	|	2983.9556	|	3431.838	|	2040	|	4584	|	FALSE |	normal	|
|	2021-01-14 	|	rest_of_week 	|	3059	|	    TRUE 	|	3226.592	|	3005.9247	|	3463.459	|	2082	|	4607	|	FALSE |	normal	|
|	2021-01-15 	|	rest_of_week 	|	3394	|	    TRUE 	|	3253.333	|	3028.0111	|	3495.421	|	2078	|	4761	|	FALSE |	normal	|
|	2021-01-16 	|	     weekend 	|	3094	|	    TRUE 	|	2190.757	|	1997.2333	|	2403.033	|	1401	|	3103	|	FALSE |	normal	|
|	2021-01-17 	|	     weekend 	|	1923	|	    TRUE 	|	2208.914	|	2012.2072	|	2424.850	|	1397	|	3187	|	FALSE |	normal	|
|	2021-01-18 	|	      monday 	|	2166	|	    TRUE 	|	2139.759	|	1897.3258	|	2413.170	|	1368	|	3178	|	FALSE |	normal	|
|	2021-01-19 	|	rest_of_week 	|	3819	|	    TRUE 	|	3362.531	|	3117.5679	|	3626.743	|	2168	|	4810	|	FALSE |	normal	|
|	2021-01-20 	|	rest_of_week 	|	4123	|	    TRUE 	|	3390.399	|	3140.2694	|	3660.452	|	2209	|	4928	|	FALSE |	normal	|
|	2021-01-21 	|	rest_of_week 	|	4151	|	    TRUE 	|	3418.497	|	3163.0995	|	3694.517	|	2214	|	4821	|	FALSE |	normal	|
|	2021-01-22 	|	rest_of_week 	|	3847	|	    TRUE 	|	3446.829	|	3186.0600	|	3728.941	|	2223	|	5015	|	FALSE |	normal	|
|	2021-01-23 	|	     weekend 	|	3593	|	    TRUE 	|	2321.055	|	2103.8608	|	2560.672	|	1488	|	3375	|	 TRUE |	increase	|
|	2021-01-24 	|	     weekend 	|	2269	|	    TRUE 	|	2340.291	|	2119.4468	|	2584.148	|	1442	|	3347	|	FALSE	|	normal	|
|	2021-01-25 	|	      monday 	|	2231	|	    TRUE 	|	2267.024	|	2001.0660	|	2568.330	|	1423	|	3313	|	FALSE	|	normal	|
|	2021-01-26 	|	rest_of_week 	|	4413	|	    TRUE 	|	3562.522	|	3279.2427	|	3870.272	|	2248	|	5189	|	FALSE	|	normal	|
|	2021-01-27 	|	rest_of_week 	|	3789	|	    TRUE 	|	3592.047	|	3302.8822	|	3906.528	|	2253	|	5123	|	FALSE	|	normal	|
|	2021-01-28 	|	rest_of_week 	|	4159	|	    TRUE 	|	3621.817	|	3326.6629	|	3943.158	|	2241	|	5154	|	FALSE	|	normal	|
|	2021-01-29 	|	rest_of_week 	|	3724	|	    TRUE 	|	3651.833	|	3350.5862	|	3980.165	|	2320	|	5402	|	FALSE	|	normal	|


![unnamed-chunk-6-1](https://user-images.githubusercontent.com/81782228/126196136-001668fa-0634-4c98-b8f1-0ee9b411bc38.png)



 
This result and method was forked from https://github.com/reconhub/trendbreaker and modified to apply PAHO regional data. 

