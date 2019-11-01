# ddfeff
 Directional Distance Function for Efficiency/Productivity Analysis in Stata
 
 ## Installing within Stata
 
```
 net install sbmeff, from("https://raw.githubusercontent.com/kerrydu/ddfeff/master/")
```

 Alternatively, download the zipfile and unzip it to your computer disk. 
```
   copy https://codeload.github.com/kerrydu/ddfeff/zip/master ddfeff-master.zip
   unzipfile ddfeff-master.zip
   net install ddfeff, from(`c(pwd)'/ddfeff-master)
```