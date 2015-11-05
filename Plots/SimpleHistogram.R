
input_table<-WorldPhones
input_region<-"Europe"

{{
# Render a barplot
    barplot(input_table[,input_region]*1000, 
            main=input_region,
            ylab="Number of Telephones",
            xlab="Year")
}}
