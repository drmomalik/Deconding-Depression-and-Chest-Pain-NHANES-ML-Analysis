## Method for extracting if respondents are on anti-depressants 

test <- nhanes("RXQ_RX_E", translated = F)
ref <- nhanes("RXQ_DRUG", translated = F)
# Filter dataframe to include rows where the column (e.g., RXDRUG) contains "antidepressant"
# Filter the dataframe where "antidepressant" appears in any of the specified columns
ref <- ref[
  grepl("antidepressant", ref$RXDDCN1B, ignore.case = TRUE) |
    grepl("antidepressant", ref$RXDDCN1C, ignore.case = TRUE) |
    grepl("antidepressant", ref$RXDDCN2A, ignore.case = TRUE) |
    grepl("antidepressant", ref$RXDDCN2B, ignore.case = TRUE) |
    grepl("antidepressant", ref$RXDDCN2C, ignore.case = TRUE),
]
# create a vector with all antidepressant ID 
id <- ref$RXDDRGID

# List of values to check
codes <- c("a70010", "c00249", "c00250", "d00144", "d00145", "d00146", "d00181", "d00217", 
           "d00236", "d00259", "d00395", "d00873", "d00874", "d00875", "d00876", "d00877", 
           "d00880", "d00882", "d00883", "d00884", "d00976", "d03157", "d03181", "d03804", 
           "d03808", "d04025", "d04332", "d04812", "d05355", "d06635", "d07113", "d07740", 
           "d08114", "d08125", "h00035")

# Create MEDDEP column
test$MEDDEP <- ifelse(test$RXDDRGID %in% codes, 1, 0)

# select only SEQN and MEDDEP 
test <- test %>% select(SEQN, MEDDEP)

#collapse over unqiue SEQN 
test <- test %>%
  group_by(SEQN) %>%               # Group by SEQN
  summarize(MEDDEP = max(MEDDEP)) # Take the maximum value of MEDDEP for each SEQN


## Only works for later cyceles (alt method)
##

test <- nhanes("P_RXQ_RX", translated = F)
test$MEDDEP <- ifelse(test$RXDRSC1 == "F32.9" |test$RXDRSC1 == "F33.9"|test$RXDRSC1 == "F39" |
                        test$RXDRSC2 == "F32.9" |test$RXDRSC2 == "F33.9"|test$RXDRSC2 == "F39" |
                        test$RXDRSC3 == "F32.9" |test$RXDRSC3 == "F33.9"|test$RXDRSC3 == "F39", 1, 0)
test2 <- test %>% select(SEQN, MEDDEP)

#collapse over unqiue SEQN 
test3 <- test2 %>%
  group_by(SEQN) %>%               # Group by SEQN
  summarize(MEDDEP = max(MEDDEP)) # Take the maximum value of MEDDEP for each SEQN