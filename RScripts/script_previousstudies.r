# previous measurements of MW mass

# comparisons (r, M, lower, upper)
Kochanek1996 = c(50, 0.51, 0.40, 0.64)
Kochanek1996woLeoI = c(50, 0.39, 0.32, 0.55)
WilkinsonEvans = c(50, 0.54, 0.18, 0.56)
# Battaglia et al 2005 
Sakamoto2003 = c(50, 0.55, 0.52, 0.551)
Battaglia2005 = c(120, 0.54, 0.40, 0.74)
Xue2008 = c(60, 0.4, 0.33, 0.47)
Gnedin2010 = c(80, 0.69, 0.57, 0.99)
Watkins2010 = c(100, 1.38, 0.99, 1.77)
Watkins2010woDraco = c(100, 0.69, 0.5, 0.88)
McMillan2011 = c(60, 0.59, 0.54, 0.64)
# Deason2012a = c(50, 0.385, 0.33, 0.44)
Deason2012a = c(50, 0.33, 0.29, 0.37)
Deason2012b = c(150, 0.75, 0.5, 1)
Kafle2012 = c(25, 0.21, 0.18, 0.25)
Gibbons2014 = c(50, 0.29, 0.24, 0.34)
Gibbons2014100=c(100, 0.40, 0.33, 0.47)
Kupper2015 = c(19, 0.21, 0.17, 0.25)
EHW2015 = c(125, 1.203, 1.121, 1.328)
McMillan2017 = c(100, 0.82, 0.71, 0.93)
Posti2019 = c(20, 0.191, 0.176, 0.208)
Sohn2018b = c(39.5, 0.61, 0.49, 0.79)
MalhanIbata2018 = c(14.5, 0.175, 0.17, 0.181)
Vasiliev2019 = c(50, 0.54, 0.46, 0.65)
Vasiliev2019100 = c(100, 0.85, 0.65, 1.18)
Watkins2019 = c(21.1, 0.22, 0.19, 0.26)
Watkins2019HST = c(39.5, 0.44, 0.38, 0.51)

# make plot for Banting application
namecompare = c("Kochanek 1996",
                "Kochanek 1996 (w/0 Leo I)",
                "Wilkinson & Evans 1999",
                "Sakamoto et al 2003",
                "Battaglia et al 2005",
                "Xue et al 2008", 
                "Gnedin et al 2010",
                "Watkins et al 2010",
                "Watkins et al 2010 (without Draco)",
                "McMillan 2011",
                "Deason et al 2012a",
                "Deason et al 2012b",
                "Kafle et al 2012",
                "Gibbons et al 2014",
                "Gibbons et al 2014",
                "EHW 2015 (with DGs)",
                "Kupper et al 2015",
                "McMillan 2017",
                "Malhan & Ibata 2018 (GD 1, Gaia)",
                "Sohn et al 2018b",
                "Posti & Helmi 2019 (Gaia)",
                "Vasiliev 2019 (Gaia)",
                "Vasiliev 2019 (at 100kpc)",
                "Watkins et al 2019 (Gaia)",
                "Watkins et al 2019 (Gaia + HST)")

ncompare = length(namecompare)

comparethese = matrix(data = c(Kochanek1996, Kochanek1996woLeoI, WilkinsonEvans, Sakamoto2003, Battaglia2005, Xue2008, Gnedin2010, 
                               Watkins2010, Watkins2010woDraco, McMillan2011, Deason2012a, Deason2012b, Kafle2012, Gibbons2014, Gibbons2014100, EHW2015, Kupper2015, McMillan2017, MalhanIbata2018,  Sohn2018b, Posti2019, Vasiliev2019, Vasiliev2019100, Watkins2019, Watkins2019HST), nrow = ncompare, byrow = TRUE)

pchcompare= seq(1:ncompare)

colcompare = rainbow(n = ncompare, start=0.05, end=1, alpha=0.85)

previousstudies = as.data.frame(comparethese) 
colnames(previousstudies) = c("r", "Mr", "lower", "upper" )

previousstudies$Study = namecompare
