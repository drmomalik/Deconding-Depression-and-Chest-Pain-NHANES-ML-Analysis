library(flowchart)

df <- subsetd

fc1 <- df %>% 
  as_fc(label = "NHANES Cardiovascular Health Survey participants between 2005-2020 (pre-pandemic)
        (All participants 40 years of age or older)", text_fs = 11, text_padding = 1) %>% 
  fc_filter(N = 7680, text_padding = 1, text_fs = 11, text_fs_exc = 11, offset_exc = 0.15, text_padding_exc = 1,
            label="Total respondents who have ever experienced chest pain", 
            text_pattern_exc = "Responded 'No' to experiencing chest pain: 21123 (72.8%)
            Refused to respond: 3 (<0.001%)
            Missing response: 20 (<0.001%)"
            , show_exc = TRUE, just_exc = "left", text_fface_exc = 1) %>% 
  fc_split(N = c(1118, 5672, 890), text_padding = 1.5, text_fs = 11,
            label = c("Total respondents with history of depression", "Total respondents with no history of depression", 
                      "No response or missing data for depression variable"), 
            ) %>% 
  fc_split(N = c(301, 817, 1429, 4243), sel_group=c("group 1","group 2"),
  label = c("History of Coronary Artery Disease/MI/Angina", "No History of Coronary Artery Disease/MI/Angina"
  ),
           text_fs = 10, text_padding = 0.5) %>% 
  fc_draw(title = "NHANES Chest-Pain and Depression Survey Population", arrow_type = "open", title_y = 0.95)




fc2 <- as_fc(N = 7000, label = "each respondent localized pain in up to 8 different body regions",
             text_pattern = "For n = 1531, {label}",
             text_fs = 12, text_padding = 1.5) %>% 
  fc_filter(N = 2487, label = "Total observations of chest pain by location", text_padding = 1.5,
            text_color_exc = "red", border_color_exc = "red", direction_ex = "right", 
            show_exc = TRUE, just = "centre", just_exc = "left", text_pattern = "{label} 
            n* = {n}",
            text_pattern_exc = "{label_exc}",
            offset_exc = 0.15, text_padding_exc = 1, text_fs_exc = 10, text_fs = 12,
            label_exc = "# of respondents with # of body regions selected:
                  n = 916 with 1 body region,
                  n = 381 with 2 body regions,
                  n = 164 with 3 body regions,
                  n = 51 with 4 body regions,
                  n = 10 with 5 body regions, 
                  n = 4 with 6 body regions, 
                  n = 1 with 7 body regions,
                  n = 4 with 8 body regions", border_color = "red",text_color = "red"
  ) %>% 
  fc_split(N = c(124, 285, 153, 883, 299, 572, 126, 45), 
         label = c("Right Arm","Right Chest",
                   "Neck","Upper Sternum","Lower Sternum",
                   "Left Chest", "Left arm", "Epigastric"), text_fs = 12
      , text_padding = 0.5) %>% 
  fc_draw(title = "Total Observations of Chest Pain by location in NHANES sample",
          arrow_type = "open", title_y = 0.95)


  
