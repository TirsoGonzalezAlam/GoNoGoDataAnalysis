# -*- coding: utf-8 -*-
"""
Created on Mon Mar 07 16:41:08 2016

@author: trga500
EFFICIENCY SCORES: RT DIVIDED BY ACCURACY (BY A FRACTION)
Absolute paths are currently not working. You need to copy and paste your \
local results directory in line 27
ADD ERROR BARS TO ACCURACY PLOT
"""

#Import things we need
import os, csv
import numpy as np
import scipy
from scipy import stats
from matplotlib import pyplot

#Implement absolute path and change directory to the script's location
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

#Load the directory to analyse
parent = os.listdir(dname)

#Prepare csv file for output

f = open(r'C:\Users\trga500\Desktop\Data\Collection\Week 6\GoNoGo\Results\WeeklyAnalysis.csv', 'wb+')
filewriter = csv.writer(f, delimiter=',')
f.write('Part_ID, Word_RT_Mean, Word_SD, Word_Acc, Word_Eff, Pic_RT_Mean, \
Pic_SD, Pic_Acc, Pic_Eff, Easy_RT_Mean, Easy_SD, Easy_Acc, Easy_Eff, \
Hard_RT_Mean, Hard_SD, Hard_Acc, Hard_Eff\n')


#Initialise empty variables and lists we need

g_word_rt = []
g_pic_rt = []
g_easy_rt = []
g_hard_rt = []
g_word_acc = []
g_pic_acc = []
g_easy_acc = []
g_hard_acc = []
g_word_eff = []
g_pic_eff = []
g_easy_eff = []
g_hard_eff = []
g_word_corrperc = []
g_pic_corrperc = []
g_easy_corrperc = []
g_hard_corrperc = []
w_mean = []
p_mean = []
e_mean = []
h_mean = []
aw_mean = []
ap_mean = []
ae_mean = []
ah_mean = []

g_word_counter = 0
g_pic_counter = 0
g_easy_counter = 0
g_hard_counter = 0

#Working version of printing a header in a very dumb way
header = '| Participant |        Word (s)        |        Image (s)       |   \
     Easy (s)        |        Hard (s)       |\n|     ID      | Mean | Stddev \
| % Corr | Mean | Stddev | % Corr | Mean | Stddev | % Corr | Mean | \
Stddev | % Corr|'
g_header = '|            Word (g)             |            Image (g)            \
|             Easy (g)            |              Hard (g)          |\n| Mean \
| Stddev | % Corr | Eff_Sc | Mean | Stddev | % Corr | Eff_Sc | Mean | Stddev \
| % Corr | Eff_Sc | Mean | Stddev | % Corr| Eff_Sc |'
print '=======================================================================\
==========================================='
print '|                                            Individual Statistics     \
                                          |'
print '=======================================================================\
==========================================='
print header



#Now we fill in the rows
for i in parent:
    word_rt = []
    pic_rt = []
    easy_rt = []
    hard_rt = []
    word_acc = []
    pic_acc = []
    easy_acc = []
    hard_acc = []
    word_counter = 0
    pic_counter = 0
    easy_counter = 0
    hard_counter = 0
    if i.endswith('.csv'):
        g = open(i)
        for line in g:
            #Prepare the lines
            line2 = line.strip()
            data = line2.split(',')
            if len(data) < 3:
                continue
            else:
                #Load the RT and ACC data into the appropriate list                
                if data[2] == 'words easy' or data[2] == 'words hard':
                    #for each trial
                    if data[3] == 'NO GO':
                        word_counter = word_counter + 1
                        word_acc.append(data[6])
                    if data[3] == 'GO':
                        if data[8] != 'nan':
                            word_rt.append(float(data[8]))
                            
                    #and for the global statistics
                    if data[3] == 'NO GO':
                        g_word_counter = g_word_counter + 1
                        g_word_acc.append(data[6])
                    if data[3] == 'GO':
                        if data[8] != 'nan':
                            g_word_rt.append(float(data[8]))
                                        
                elif data[2] == 'image easy' or data[2] == 'image hard':
                    #for each trial
                    if data[3] == 'NO GO':
                        pic_counter = pic_counter + 1
                        pic_acc.append(data[6])
                    if data[3] == 'GO':
                        if data[8] != 'nan':
                            pic_rt.append(float(data[8]))
                    
                    #and for the global statistics
                    if data[3] == 'NO GO':
                        g_pic_counter = g_pic_counter + 1
                        g_pic_acc.append(data[6])
                    if data[3] == 'GO':
                        if data[8] != 'nan':
                            g_pic_rt.append(float(data[8]))
                   
                elif data[2] == 'scrambled words easy' or data[2] == \
                'scrambled pics easy':
                    #for each trial
                    if data[3] == 'NO GO':
                        easy_counter = easy_counter + 1
                        easy_acc.append(data[6])
                    if data[3] == 'GO':
                        if data[8] != 'nan':
                            easy_rt.append(float(data[8]))
                    
                    #and for the global statistics
                    if data[3] == 'NO GO':
                        g_easy_counter = g_easy_counter + 1
                        g_easy_acc.append(data[6])
                    if data[3] == 'GO':
                        if data[8] != 'nan':
                            g_easy_rt.append(float(data[8]))
                    
                    
                elif data[2] == 'scrambled words hard' or data[2] == \
                'scrambled pics hard':
                    #for each trial
                    if data[3] == 'NO GO':
                        hard_counter = hard_counter + 1
                        hard_acc.append(data[6])
                    if data[3] == 'GO':
                        if data[8] != 'nan':
                            hard_rt.append(float(data[8]))
                    
                    #and for the global statistics
                    if data[3] == 'NO GO':
                        g_hard_counter = g_hard_counter + 1
                        g_hard_acc.append(data[6])
                    if data[3] == 'GO':
                        if data[8] != 'nan':
                            g_hard_rt.append(float(data[8]))
        g.close()
                    
                    
        #Once we have this data, we can calculate individual results
        #Calculate mean and %correct for individual participant
        #Calculate ID
        part_id = i
        
        #Word
        x = stats.itemfreq(word_acc)
        try:
            corr_perc_w = (float(x[1][1])/word_counter)*100
        except IndexError:
            corr_perc_w = 100
        mean_rt_w = np.mean(word_rt)
        mean_sd_w = np.std(word_rt)
        w_mean.append(mean_rt_w)
        aw_mean.append(corr_perc_w)
        word_eff = mean_rt_w/(corr_perc_w/100)
        g_word_eff.append(word_eff)
        g_word_corrperc.append(corr_perc_w)
        
        #Picture
        y = stats.itemfreq(pic_acc)
        try:
            corr_perc_p = (float(y[1][1])/pic_counter)*100
        except IndexError:
            corr_perc_p = 100
        mean_rt_p = np.mean(pic_rt)
        mean_sd_p = np.std(pic_rt)
        p_mean.append(mean_rt_p)
        ap_mean.append(corr_perc_p)
        pic_eff = mean_rt_p/(corr_perc_p/100)
        g_pic_eff.append(pic_eff)
        g_pic_corrperc.append(corr_perc_p)
        
        #Easy
        z = stats.itemfreq(easy_acc)
        try:
            corr_perc_e = (float(z[1][1])/easy_counter)*100
        except IndexError:
            corr_perc_e = 100
        mean_rt_e = np.mean(easy_rt)
        mean_sd_e = np.std(easy_rt)
        e_mean.append(mean_rt_e)
        ae_mean.append(corr_perc_e)
        easy_eff = mean_rt_e/(corr_perc_e/100)
        g_easy_eff.append(easy_eff)
        g_easy_corrperc.append(corr_perc_e)
        
        #Hard
        a = stats.itemfreq(hard_acc)
        try:
            corr_perc_h = (float(a[1][1])/hard_counter)*100
        except IndexError:
            corr_perc_h = 100
        mean_rt_h = np.mean(hard_rt)
        mean_sd_h = np.std(hard_rt)
        h_mean.append(mean_rt_h)
        ah_mean.append(corr_perc_h)
        hard_eff = mean_rt_h/(corr_perc_h/100)
        g_hard_eff.append(hard_eff)
        g_hard_corrperc.append(corr_perc_h)
        
        #Write results into a csv file
        f.write('%s, %.4f, %.2f, %.2f, %.2f, %.4f, %.2f, %.2f, %.2f, %.4f, \
        %.2f, %.2f, %.2f, %.4f, %.2f, %.2f, %.2f\n'%(part_id[:2].rstrip('_'), \
        mean_rt_w, mean_sd_w, corr_perc_w, word_eff, mean_rt_p, mean_sd_p, \
        corr_perc_p, pic_eff, mean_rt_e, mean_sd_e, corr_perc_e, easy_eff, \
        mean_rt_h, mean_sd_h, corr_perc_h, hard_eff))
        
        #Print each row (Currently not printing efficiency scores)
        if corr_perc_w == 100.0:
            print '|      %s     | %.2f |  %.2f  | %.1f  | %.2f |  %.2f  |  \
%.2f | %.2f |  %.2f  | %.2f  | %.2f |  %.2f  |  %.2f|'%(part_id[:2], \
            mean_rt_w, mean_sd_w, corr_perc_w, mean_rt_p, mean_sd_p, \
            corr_perc_p, mean_rt_e, mean_sd_e, corr_perc_e, mean_rt_h, \
            mean_sd_h, corr_perc_h)
        elif corr_perc_e != 100.0:
            print '|      %s     | %.2f |  %.2f  | %.2f  | %.2f |  %.2f  |  \
%.2f | %.2f |  %.2f  | %.2f  | %.2f |  %.2f  |  %.2f|'%(part_id[:2], \
            mean_rt_w, mean_sd_w, corr_perc_w, mean_rt_p, mean_sd_p, \
            corr_perc_p, mean_rt_e, mean_sd_e, corr_perc_e, mean_rt_h, \
            mean_sd_h, corr_perc_h)
        else:
            print '|      %s     | %.2f |  %.2f  | %.2f  | %.2f |  %.2f  |  \
%.2f | %.2f |  %.2f  | %.1f  | %.2f |  %.2f  |  %.2f|'%(part_id[:2], \
            mean_rt_w, mean_sd_w, corr_perc_w, mean_rt_p, mean_sd_p, \
            corr_perc_p, mean_rt_e, mean_sd_e, corr_perc_e, mean_rt_h, \
            mean_sd_h, corr_perc_h)        

print 'Note. The RT data comes from the Go trials, and the accuracy from the \
No Go trials'

#Calculate z scores per participant
#RT
z_score_rt_w = scipy.stats.zscore(w_mean)
z_score_rt_p = scipy.stats.zscore(p_mean)
z_score_rt_e = scipy.stats.zscore(e_mean)
z_score_rt_h = scipy.stats.zscore(h_mean)
#ACC
z_score_rt_aw = scipy.stats.zscore(aw_mean)
z_score_rt_ap = scipy.stats.zscore(ap_mean)
z_score_rt_ae = scipy.stats.zscore(ae_mean)
z_score_rt_ah = scipy.stats.zscore(ah_mean)
#EFF
z_score_rt_ew = scipy.stats.zscore(g_word_eff)
z_score_rt_ep = scipy.stats.zscore(g_pic_eff)
z_score_rt_ee = scipy.stats.zscore(g_easy_eff)
z_score_rt_eh = scipy.stats.zscore(g_hard_eff)

#Write the z scores results in the results file
f.write('\n')
f.write('z Scores \n')
f.write('ID, z_RT_word, z_RT_pic, z_RT_easy, z_RT_hard, z_ACC_word, z_ACC_pic, \
z_ACC_easy, z_ACC_hard, z_EFF_word, z_EFF_pic, z_EFF_easy, z_EFF_hard \n')
for i in range(0, len(z_score_rt_w)):
    f.write(('%s, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, \
    %.3f, %.3f \n')%(int(i)+1, z_score_rt_w[i], z_score_rt_p[i], \
    z_score_rt_e[i], z_score_rt_h[i], z_score_rt_aw[i], z_score_rt_ap[i], \
    z_score_rt_ae[i], z_score_rt_ah[i], z_score_rt_ew[i], z_score_rt_ep[i], \
    z_score_rt_ee[i], z_score_rt_eh[i]))

#Once we're done, we can calculate the group statistics
#Word
g_x = stats.itemfreq(g_word_acc)
g_corr_perc_w = (float(g_x[1][1])/g_word_counter)*100
g_mean_rt_w = np.mean(g_word_rt)
g_mean_sd_w = np.std(g_word_rt)
g_mean_corrperc_word = np.mean(g_word_corrperc)
g_mean_eff_w = np.mean(g_word_eff)
normality_test_w_RT = stats.mstats.normaltest(w_mean)
normality_test_w_ACC = stats.mstats.normaltest(aw_mean)
normality_test_w_EFF = stats.mstats.normaltest(g_word_eff)
StdErr_RT_w = stats.sem(g_word_rt)
StdErr_CP_w = stats.sem(g_word_corrperc)
StdErr_EFF_w = stats.sem(g_word_eff)

 #Picture
g_y = stats.itemfreq(g_pic_acc)
g_corr_perc_p = (float(g_y[1][1])/g_pic_counter)*100
g_mean_rt_p = np.mean(g_pic_rt)
g_mean_sd_p = np.std(g_pic_rt)
g_mean_corrperc_pic = np.mean(g_pic_corrperc)
g_mean_eff_p = np.mean(g_pic_eff)
normality_test_p_RT = stats.mstats.normaltest(p_mean)
normality_test_p_ACC = stats.mstats.normaltest(ap_mean)
normality_test_p_EFF = stats.mstats.normaltest(g_pic_eff)
StdErr_RT_p = stats.sem(g_pic_rt)
StdErr_CP_p = stats.sem(g_pic_corrperc)
StdErr_EFF_p = stats.sem(g_pic_eff)

 #Easy
g_z = stats.itemfreq(g_easy_acc)
g_corr_perc_e = (float(g_z[1][1])/g_easy_counter)*100
g_mean_rt_e = np.mean(g_easy_rt)
g_mean_sd_e = np.std(g_easy_rt)
g_mean_corrperc_easy = np.mean(g_easy_corrperc)
g_mean_eff_e = np.mean(g_easy_eff)
normality_test_e_RT = stats.mstats.normaltest(e_mean)
normality_test_e_ACC = stats.mstats.normaltest(ae_mean)
normality_test_e_EFF = stats.mstats.normaltest(g_easy_eff)
StdErr_RT_e = stats.sem(g_easy_rt)
StdErr_CP_e = stats.sem(g_easy_corrperc)
StdErr_EFF_e = stats.sem(g_easy_eff)

 #Hard
g_a = stats.itemfreq(g_hard_acc)
g_corr_perc_h = (float(g_a[1][1])/g_hard_counter)*100
g_mean_rt_h = np.mean(g_hard_rt)
g_mean_sd_h = np.std(g_hard_rt)
g_mean_corrperc_hard = np.mean(g_hard_corrperc)
g_mean_eff_h = np.mean(g_hard_eff)
normality_test_h_RT = stats.mstats.normaltest(h_mean)
normality_test_h_ACC = stats.mstats.normaltest(ah_mean)
normality_test_h_EFF = stats.mstats.normaltest(g_hard_eff)
StdErr_RT_h = stats.sem(g_hard_rt)
StdErr_CP_h = stats.sem(g_hard_corrperc)
StdErr_EFF_h = stats.sem(g_hard_eff)

#And now, the inferential statistics
anova_rt = scipy.stats.f_oneway(w_mean, p_mean, e_mean, h_mean)
anova_acc = scipy.stats.f_oneway(aw_mean, ap_mean, ae_mean, ah_mean)
anova_eff = scipy.stats.f_oneway(g_word_eff, g_pic_eff, g_easy_eff, g_hard_eff)

#Print group table
print ''
print ''
print '==================================================================\
======================================================================'
print '|                                                            Group \
Statistics                                                          |'
print '====================================================================\
===================================================================='
print g_header
print '| %.2f |  %.2f  | %.2f  |  %.2f  | %.2f |  %.2f  |  %.2f |  %.2f  | \
%.2f |  %.2f  |  %.2f |  %.2f  | %.2f |  %.2f  |  %.2f|  %.2f  | \
'%(g_mean_rt_w, g_mean_sd_w, g_corr_perc_w, g_mean_eff_w, g_mean_rt_p, \
g_mean_sd_p, g_corr_perc_p, g_mean_eff_p, g_mean_rt_e, g_mean_sd_e, \
g_corr_perc_e, g_mean_eff_e, g_mean_rt_h, g_mean_sd_h, g_corr_perc_h, \
g_mean_eff_h)
print ''
print ''
print 'The following figure shows the distribution of RT means per condition: \
\n (blue = word, red = picture, yellow = easy box, green = hard box)'

f.write('\nGROUP RESULTS\n')
f.write('Part_ID, Word_RT_Mean, Word_SD, Word_Acc, Word_Eff, Pic_RT_Mean, \
Pic_SD, Pic_Acc, Pic_Eff, Easy_RT_Mean, Easy_SD, Easy_Acc, Easy_Eff, \
Hard_RT_Mean, Hard_SD, Hard_Acc, Hard_Eff\n')
f.write('GROUP, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f'%(g_mean_rt_w, g_mean_sd_w, g_corr_perc_w, g_mean_eff_w, g_mean_rt_p, g_mean_sd_p, g_corr_perc_p, g_mean_eff_p, g_mean_rt_e, g_mean_sd_e, g_corr_perc_e, g_mean_eff_e, g_mean_rt_h, g_mean_sd_h, g_corr_perc_h, g_mean_eff_h))

#==============================================================================
# DISTRIBUTION PLOTS
#==============================================================================

#RT plot
pyplot.figure(1,figsize=(8,6),dpi=300)
pyplot.subplot(221)
pyplot.ylabel('Frequency')
pyplot.xlabel('Reaction time')
pyplot.axis([0.3, 0.8, 0, 6])
pyplot.hist(w_mean, histtype='bar', label='Word RT mean', color='blue')
pyplot.subplot(222)
pyplot.ylabel('Frequency')
pyplot.xlabel('Reaction time')
pyplot.axis([0.3, 0.8, 0, 6])
pyplot.hist(p_mean, histtype='bar', label='Picture RT mean', color='red')
pyplot.subplot(223)
pyplot.ylabel('Frequency')
pyplot.xlabel('Reaction time')
pyplot.axis([0.3, 0.8, 0, 6])
pyplot.hist(e_mean, histtype='bar', label='Easy RT mean', color='yellow')
pyplot.subplot(224)
pyplot.ylabel('Frequency')
pyplot.xlabel('Reaction time')
pyplot.axis([0.3, 0.8, 0, 6])
pyplot.hist(h_mean, histtype='bar', label='Hard RT mean', color='green')
pyplot.show()

#Accuracy plot
print ''
print ''
print 'The following figure shows the distribution of Accuracy means per \
condition: \n (blue = word, red = picture, yellow = easy box, green = hard box)'
pyplot.figure(2,figsize=(8,6),dpi=300)
pyplot.subplot(221)
pyplot.ylabel('Frequency')
pyplot.xlabel('Accuracy')
pyplot.axis([20, 100, 0, 7])
pyplot.hist(aw_mean, histtype='bar', label='Word ACC mean', color='blue')
pyplot.subplot(222)
pyplot.ylabel('Frequency')
pyplot.xlabel('Accuracy')
pyplot.axis([20, 100, 0, 7])
pyplot.hist(ap_mean, histtype='bar', label='Picture ACC mean', color='red')
pyplot.subplot(223)
pyplot.ylabel('Frequency')
pyplot.xlabel('Accuracy')
pyplot.axis([20, 100, 0, 7])
pyplot.hist(ae_mean, histtype='bar', label='Easy ACC mean', color='yellow')
pyplot.subplot(224)
pyplot.ylabel('Frequency')
pyplot.xlabel('Accuracy')
pyplot.axis([20, 100, 0, 7])
pyplot.hist(ah_mean, histtype='bar', label='Hard ACC mean', color='green')
pyplot.show()

#Efficiency Plot
print ''
print ''
print 'The following figure shows the distribution of Efficiency Scores per\
 condition: \n (blue = word, red = picture, yellow = easy box, green = hard box)'
pyplot.figure(3,figsize=(8,6),dpi=300)
pyplot.subplot(221)
pyplot.ylabel('Frequency')
pyplot.xlabel('Efficiency')
pyplot.axis([0.4, 2.0, 0, 5])
pyplot.hist(g_word_eff, histtype='bar', label='Word EFF mean', color='blue')
pyplot.subplot(222)
pyplot.ylabel('Frequency')
pyplot.xlabel('Efficiency')
pyplot.axis([0.4, 2.0, 0, 5])
pyplot.hist(g_pic_eff, histtype='bar', label='Picture EFF mean', color='red')
pyplot.subplot(223)
pyplot.ylabel('Frequency')
pyplot.xlabel('Efficiency')
pyplot.axis([0.4, 2.0, 0, 8])
pyplot.hist(g_easy_eff, histtype='bar', label='Easy EFF mean', color='yellow')
pyplot.subplot(224)
pyplot.ylabel('Frequency')
pyplot.xlabel('Efficiency')
pyplot.axis([0.4, 2.0, 0, 12])
pyplot.hist(g_hard_eff, histtype='bar', label='Hard EFF mean', color='green')
pyplot.show()

#==============================================================================
# GROUP MEANS PLOT
#==============================================================================

###RT###
# Get means we need to plot
RT_Means = g_mean_rt_w, g_mean_rt_p, g_mean_rt_e, g_mean_rt_h
StdErr_RT = [StdErr_RT_w, StdErr_RT_p, StdErr_RT_e, StdErr_RT_h]
width = 0.7
ind = [1, 2, 3, 4]
labels = ["Word", "Picture", "Easy", "Hard"]
# Plot them
f4 = pyplot.figure(4,figsize=(8,6),dpi=500)
fig, ax = pyplot.subplots()
f4 = pyplot.ylabel('Reaction Time (seconds)')
f4 = pyplot.xlabel('Condition')
f4 = ax.set_ylim([0.46, 0.54])
f4 = ax.set_xlim([0.5, 4.5])
f4 = pyplot.xticks(ind, labels, rotation='horizontal')
f4 = pyplot.tick_params(axis='both', which='both', bottom='on', top='off',\
                   left='on', right='off')
f4 = pyplot.bar([1, 2, 3, 4], RT_Means, width, color='w', yerr=StdErr_RT, \
align='center')
f4 = pyplot.title('Reaction time group means per condition')
f4 = pyplot.savefig(r'Results\RT.png')
f4 = pyplot.show()

###ACC###
# Get means we need to plot
ACC_Means = g_corr_perc_w, g_corr_perc_p, g_corr_perc_e, g_corr_perc_h
StdErr_ACC = [StdErr_CP_w, StdErr_CP_p, StdErr_CP_e, StdErr_CP_h]
# Plot them
f5 = pyplot.figure(5,figsize=(8,6),dpi=500)
fig, ax = pyplot.subplots()
f5 = pyplot.ylabel('Accuracy (%)')
f5 = pyplot.xlabel('Condition')
f5 = ax.set_ylim([60.0, 100.0])
f5 = ax.set_xlim([0.5, 4.5])
f5 = pyplot.xticks(ind, labels, rotation='horizontal')
f5 = pyplot.tick_params(axis='both', which='both', bottom='on', top='off',\
                   left='on', right='off')
f5 = pyplot.bar([1, 2, 3, 4], ACC_Means, width, color='w', align='center', \
yerr=StdErr_ACC)
f5 = pyplot.title('Accuracy group means per condition')
f5 = pyplot.savefig(r'Results\ACC.png')
f5 = pyplot.show()

###EFFICIENCY###
# Get means we need to plot
EFF_Means = g_mean_eff_w, g_mean_eff_p, g_mean_eff_e, g_mean_eff_h
StdErr_EFF = [StdErr_EFF_w, StdErr_EFF_p, StdErr_EFF_e, StdErr_EFF_h]

# Plot them
f6 = pyplot.figure(6,figsize=(8,6),dpi=500)
fig, ax = pyplot.subplots()
f6 = pyplot.ylabel('Efficiency score')
f6 = pyplot.xlabel('Condition')
f6 = ax.set_ylim([0.5, 1.0])
f6 = ax.set_xlim([0.5, 4.5])
f6 = pyplot.xticks(ind, labels, rotation='horizontal')
f6 = pyplot.tick_params(axis='both', which='both', bottom='on', top='off',\
                   left='on', right='off')
f6 = pyplot.bar([1, 2, 3, 4], EFF_Means, width, color='w', yerr=StdErr_EFF, \
align='center')
f6 = pyplot.title('Efficiency score group means per condition')
f6 = pyplot.savefig(r'Results\EFF.png')
f6 = pyplot.show()

#pyplot.subplot(221)
#pyplot.ylabel('Frequency')
#pyplot.xlabel('Reaction time')
#pyplot.axis([0.3, 0.8, 0, 6])
#pyplot.hist(w_mean, histtype='bar', label='Word RT mean', color='blue')
#pyplot.subplot(222)
#pyplot.ylabel('Frequency')
#pyplot.xlabel('Reaction time')
#pyplot.axis([0.3, 0.8, 0, 6])
#pyplot.hist(p_mean, histtype='bar', label='Picture RT mean', color='red')
#pyplot.subplot(223)
#pyplot.ylabel('Frequency')
#pyplot.xlabel('Reaction time')
#pyplot.axis([0.3, 0.8, 0, 6])
#pyplot.hist(e_mean, histtype='bar', label='Easy RT mean', color='yellow')
#pyplot.subplot(224)
#pyplot.ylabel('Frequency')
#pyplot.xlabel('Reaction time')
#pyplot.axis([0.3, 0.8, 0, 6])
#pyplot.hist(h_mean, histtype='bar', label='Hard RT mean', color='green')
#pyplot.show()

#Print inferential statistics
print ''
print "Normality tests (D'Agostino & Pearson's test, 1973):"
print ''
print "RT"
print "There's a %.3f probability that the RT data for words are normally distributed"%(normality_test_w_RT[1])
print "There's a %.3f probability that the RT data for pictures are normally distributed"%(normality_test_p_RT[1])
print "There's a %.3f probability that the RT data for easy are normally distributed"%(normality_test_e_RT[1])
print "There's a %.3f probability that the RT data for har are normally distributed"%(normality_test_h_RT[1])
print ''
print "Accuracy"
print "There's a %.3f probability that the Accuracy data for words are normally distributed"%(normality_test_w_ACC[1])
print "There's a %.3f probability that the Accuracy data for pictures are normally distributed"%(normality_test_p_ACC[1])
print "There's a %.3f probability that the Accuracy data for easy are normally distributed"%(normality_test_e_ACC[1])
print "There's a %.3f probability that the Accuracy data for har are normally distributed"%(normality_test_h_ACC[1])
print ''
print "Efficiency"
print "There's a %.3f probability that the Efficiency data for words are normally distributed"%(normality_test_w_EFF[1])
print "There's a %.3f probability that the Efficiency data for pictures are normally distributed"%(normality_test_p_EFF[1])
print "There's a %.3f probability that the Efficiency data for easy are normally distributed"%(normality_test_e_EFF[1])
print "There's a %.3f probability that the Efficiency data for har are normally distributed"%(normality_test_h_EFF[1])
print ''
print "ANOVA"
print "The ANOVA for Accuracy returns a p value of %.3f, for Efficiency scores a p value of %.3f and for RT a p value of %.3f"%(anova_acc[1], anova_eff[1], anova_rt[1])
print "Note: The means don't come from independent data (i.e., each participant was tested in all conditions), so the assumption of independence of observations is violated. The other two assumptions (normality and homoscedasticity) seem to hold (Except for Efficiency Scores, where 2 of the distributions aren't normal). Perhaps a Repeated measures ANOVA would be more sensible"
f.close()
