K2 targets
================
J. Dotson
October 14, 2018

Goal of this document is to generate a list of unique stars observed by K2 for at least one full length campaign.

MAST has both the epic catalog and a list of all the K2 data products in their databases. The following commands were used to join these two datatables and pull the parameters of interest in the MAST CasJobs interface (<http://mastweb.stsci.edu/kplrcasjobs/SubmitJob.aspx>) The result was downloaded as a .csv file.

> SELECT INTO K2observed\_targets
> id,
> sci\_campaign,
> objtype,
> kp,
> teff,
> ep\_teff,
> em\_teff,
> logg,
> ep\_logg,
> em\_logg,
> rad,
> ep\_rad,
> em\_rad,
> mass,
> ep\_mass,
> em\_mass,
> d
> FROM k2\_science
> JOIN k2\_epic ON k2\_epic.id = k2\_science.sci\_kepler\_id
> WHERE sci\_archive\_class = 'KTL'

Start by reading in results of MAST query.

``` r
k2targets <- read_csv("data/K2observed_targets_jdotson.csv",na=c("","null"))
print(paste("Number of entries in K2targets =",nrow(k2targets)))
```

\[1\] "Number of entries in K2targets = 420949"

Remove stars observed multiple times, extended targets, and targets without teff/logg

``` r
#create list of unique targets for use in planet estimate
k2unique <- k2targets[!duplicated(k2targets$id),]
print(paste("Number of unique targets =", nrow(k2unique)))
```

\[1\] "Number of unique targets = 359265"

``` r
kable(count(k2unique,objtype),caption="unique targets by type") %>%
  kable_styling(bootstrap_options = c("striped", "condensed"),full_width=F,position="left")  
```

<table class="table table-striped table-condensed" style="width: auto !important; ">
<caption>
unique targets by type
</caption>
<thead>
<tr>
<th style="text-align:left;">
objtype
</th>
<th style="text-align:right;">
n
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
EXTENDED
</td>
<td style="text-align:right;">
34804
</td>
</tr>
<tr>
<td style="text-align:left;">
STAR
</td>
<td style="text-align:right;">
307229
</td>
</tr>
<tr>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
17232
</td>
</tr>
</tbody>
</table>
``` r
print(paste("Number of objects without Teff = ", sum(is.na(k2unique$teff))))
```

\[1\] "Number of objects without Teff = 71577"

``` r
print(paste("Number of objects without logg =", sum(is.na(k2unique$logg))))
```

\[1\] "Number of objects without logg = 71577"

Before closing the books on this one, let's just understand a little bit about the entries with out teff and logg...

``` r
need_teff <- filter(k2unique, is.na(teff))

kable(count(need_teff,objtype),caption="targets without teff/logg") %>%
  kable_styling(bootstrap_options = c("striped", "condensed"),full_width=F,position="left")
```

<table class="table table-striped table-condensed" style="width: auto !important; ">
<caption>
targets without teff/logg
</caption>
<thead>
<tr>
<th style="text-align:left;">
objtype
</th>
<th style="text-align:right;">
n
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
EXTENDED
</td>
<td style="text-align:right;">
34804
</td>
</tr>
<tr>
<td style="text-align:left;">
STAR
</td>
<td style="text-align:right;">
19541
</td>
</tr>
<tr>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
17232
</td>
</tr>
</tbody>
</table>
In general, if something doesn't have Teff and logg, then it probably isn't a star -- even if it's labeled as such. So, for the purposes of generating a list of observed stars, we can just ignore those. Save resulting list as an .RDS and .csv for use elsewhere.

``` r
k2uniquestars <- filter(k2unique, !is.na(teff))

saveRDS(k2uniquestars,"RDS/k2uniquestars.rds")
write_csv(k2uniquestars,"data/k2uniquestar.csv")
```
