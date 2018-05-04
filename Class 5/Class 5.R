#' ---
#' title: "Bioinformatics Class 5"
#' author: "Amay"


# Bioinformatics Class 5
# Plots

x <- rnorm(1000, 0)
y <- rnorm(2000, 20)

summary(x)

mean(x)
sd(x)

# Graphed as box plot
boxplot(x, y)

# visualizing histogram of x and then y 
hist(x)
hist(y)

## Section 1A from lab sheet
# read table takes a file input, everything else is default. 
# Header paramater within read.table() function takes a logical value, set to TRUE if first row contains variable names 
# check par for full detailed list of all paramaters available for modification
# pch sets shape of points, i.e. plot characters - numerous shapes assigned to numbers and vectors
# cex is area of plot point (i.e. how big point appears)
# lwd = line width
# ylim sets y axis and takes a vector input only, requires c() function
# lty is line type (dotted, dashed, solid, etc)

baby <- read.table("bggn213_05_rstats/weight_chart.txt", header = TRUE)
plot(baby, type = 'b', main = 'Baby Weights by age', xlab = 'Age (months)', ylab = 'Weight (kg)', pch = 15, cex = 1.5, lwd = 2, ylim = c(2,10), lty = 1)

##End Section 1A

## Section 1B
FC <- read.table("bggn213_05_rstats/feature_counts.txt", header = T, sep = '\t')
FC

#changing margins from default. Access default values using "par()' in console
par(mar = c(5, 11, 4, 2))

#las sets orientation of graph. 1 is normal vertical bars, 2 is horizontal bars
barplot(FC$Count, names.arg=FC$Feature, xlab = 'Number', horiz = T, las = 2)

##End Section 1B

##Section 2A
MF <- read.delim('bggn213_05_rstats/male_female_counts.txt')

#producing a vector containing 10 distinct colors
rainbow(10)

#adding colors to each bar
#barplot(MF$Count, col = c(rainbow(10)))

#customizing rainbow to specific data set by calling on number of rows
barplot(MF$Count, col = c(rainbow(nrow(MF))))

#customizing colors by sex
barplot(MF$Count, col = c(rainbow(ncol(MF))))

##End Section 2A

##Section 2B

UDX <- read.delim('bggn213_05_rstats/up_down_expression.txt')
plot(UDX$Condition1, UDX$Condition2, col = UDX$State) 

# Identifying genes by state 
# By levels
unique(UDX$State)

# number belonging to each level
table(UDX$State)

# checking palette of available colors and levels in 'State' - aligns colors to state by order in vector
palette()
levels(UDX$State)

# changing palette of available colors to desired and replotting
palette(c('blue','grey','red'))
plot(UDX$Condition1, UDX$Condition2, col = UDX$State) 

##End of Section 2B

##Section 2C

#converting methylation expression text into dataframe
Mex <- read.delim('bggn213_05_rstats/expression_methylation.txt')
plot(Mex$promoter.meth, Mex$gene.meth)

