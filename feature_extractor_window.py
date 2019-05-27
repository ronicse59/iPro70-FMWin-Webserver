import methods
import numpy as np
from fasta_reader import FASTA


def generator(input_file_name, feature_list):
    order, sequences = FASTA(input_file_name)
    print("-> Window feature set generating ...");
    newstr = []
    for s in order:
        p = sequences[s]
        p = p.upper()
        each_feature_vector = []
        # Feature 1
        # Frequencey count
        #--------------------------------------------------------------------------------------
        if(feature_list[1]):
            # on full sequence [0-80]
            a, c, g, t = methods.frequence_count(p, 'A', 'C', 'G', 'T')
            each_feature_vector.append(a)
            each_feature_vector.append(c)
            each_feature_vector.append(g)
            each_feature_vector.append(t)
            each_feature_vector.append(g+c)


            # among index[0-9]
            a,c,g,t = methods.frequence_count(p[:10], 'A', 'C', 'G', 'T')
            each_feature_vector.append(a)
            each_feature_vector.append(c)
            each_feature_vector.append(g)
            each_feature_vector.append(t)
            each_feature_vector.append(g + c)
            # among index[10-19]
            a, c, g, t = methods.frequence_count(p[10:20], 'A', 'C', 'G', 'T')
            each_feature_vector.append(a)
            each_feature_vector.append(c)
            each_feature_vector.append(g)
            each_feature_vector.append(t)
            each_feature_vector.append(g + c)
            # among index[20-34]
            a, c, g, t = methods.frequence_count(p[20:35], 'A', 'C', 'G', 'T')
            each_feature_vector.append(a)
            each_feature_vector.append(c)
            each_feature_vector.append(g)
            each_feature_vector.append(t)
            each_feature_vector.append(g + c)
            # among index[45-54]
            a, c, g, t = methods.frequence_count(p[45:55], 'A', 'C', 'G', 'T')
            each_feature_vector.append(a)
            each_feature_vector.append(c)
            each_feature_vector.append(g)
            each_feature_vector.append(t)
            each_feature_vector.append(g + c)
            # among index[55-59]
            a, c, g, t = methods.frequence_count(p[55:60], 'A', 'C', 'G', 'T')
            each_feature_vector.append(a)
            each_feature_vector.append(c)
            each_feature_vector.append(g)
            each_feature_vector.append(t)
            each_feature_vector.append(g + c)
            # among index[60-69]
            a, c, g, t = methods.frequence_count(p[60:70], 'A', 'C', 'G', 'T')
            each_feature_vector.append(a)
            each_feature_vector.append(c)
            each_feature_vector.append(g)
            each_feature_vector.append(t)
            each_feature_vector.append(g + c)
            # among index[70-80]
            a, c, g, t = methods.frequence_count(p[70:], 'A', 'C', 'G', 'T')
            each_feature_vector.append(a)
            each_feature_vector.append(c)
            each_feature_vector.append(g)
            each_feature_vector.append(t)
            each_feature_vector.append(g + c)

        # Feature 2
        # mean, variance and standard-deviation
        # --------------------------------------------------------------------------------------
        if (feature_list[2]):
            a, c, g, t = methods.frequence_count(p, 'A', 'C', 'G', 'T')
            x = [a, c, g, t]
            each_feature_vector.append(np.mean(x))
            each_feature_vector.append(np.var(x))
            each_feature_vector.append(np.std(x))

            a, c, g, t = methods.frequence_count(p[45:55], 'A', 'C', 'G', 'T')
            x = [a, c, g, t]
            each_feature_vector.append(np.mean(x))
            each_feature_vector.append(np.var(x))
            each_feature_vector.append(np.std(x))

        # Feature 3
        # G-C Skew
        # AT-GC ratio
        # --------------------------------------------------------------------------------------
        if(feature_list[3]):
            value = 0
            for s in p:
                if s == 'G':
                    value = 1
                elif s == 'C':
                    value = -1
                else:
                    value = 0
                each_feature_vector.append(value)


        # Feature 4
        # K-mar frequency count
        # K = 2,3,4,5,6
        # --------------------------------------------------------------------------------------
        # on full sequence
        if (feature_list[4]):
            each_feature_vector = each_feature_vector + methods.two_mar_frequency_count(p)
            each_feature_vector = each_feature_vector + methods.two_mar_frequency_count(p[20:35])
            each_feature_vector = each_feature_vector + methods.two_mar_frequency_count(p[45:55])
        if(feature_list[5]):
            each_feature_vector = each_feature_vector + methods.three_mar_frequency_count(p[20:35])
            each_feature_vector = each_feature_vector + methods.three_mar_frequency_count(p[45:55])
            each_feature_vector = each_feature_vector + methods.three_mar_frequency_count(p)
        if(feature_list[6]):
            each_feature_vector = each_feature_vector + methods.four_mar_frequency_count(p)
            each_feature_vector = each_feature_vector + methods.four_mar_frequency_count(p[20:35])
            each_feature_vector = each_feature_vector + methods.four_mar_frequency_count(p[45:55])
        if(feature_list[7]):
            each_feature_vector = each_feature_vector + methods.five_mar_frequency_count(p)
            each_feature_vector = each_feature_vector + methods.five_mar_frequency_count(p[20:35])
            each_feature_vector = each_feature_vector + methods.five_mar_frequency_count(p[45:55])
        if(feature_list[8]):
            each_feature_vector = each_feature_vector + methods.six_mar_frequency_count(p)
            each_feature_vector = each_feature_vector + methods.six_mar_frequency_count(p[20:35])
            each_feature_vector = each_feature_vector + methods.six_mar_frequency_count(p[45:55])


        # Feature 5
        # 2-mar and K-gap count
        # --------------------------------------------------------------------------------------
        if (feature_list[9]):
            for i in range(1, 25):
                each_feature_vector = each_feature_vector + methods.two_mar_k_gap(p, i)
            for i in range(1,12):
                each_feature_vector = each_feature_vector + methods.two_mar_k_gap(p[20:35], i)
            for i in range(1,8):
                each_feature_vector = each_feature_vector + methods.two_mar_k_gap(p[45:55], i)


        # Feature 6
        # 3-mar and right K-gap
        # --------------------------------------------------------------------------------------
        if (feature_list[10]):
            for i in range(1, 25):
                each_feature_vector = each_feature_vector + methods.three_mar_right_k_gap(p, i)
            for i in range(1, 12):
                each_feature_vector = each_feature_vector + methods.three_mar_right_k_gap(p[20:35], i)
            for i in range(1, 8):
                each_feature_vector = each_feature_vector + methods.three_mar_right_k_gap(p[45:55], i)



        # Feature 7
        # 3-mar and left K-gap
        # --------------------------------------------------------------------------------------
        if (feature_list[11]):
            for i in range(1, 25):
                each_feature_vector = each_feature_vector + methods.three_mar_left_k_gap(p, i)
            for i in range(1, 12):
                each_feature_vector = each_feature_vector + methods.three_mar_left_k_gap(p[20:35], i)
            for i in range(1, 8):
                each_feature_vector = each_feature_vector + methods.three_mar_left_k_gap(p[45:55], i)



        # Feature 8
        # Pattern matching with minmum 3 matching is acceptable
        # --------------------------------------------------------------------------------------
        if(feature_list[12]):
            st = 45
            en = 55
            tata1 = "TATAAT"
            tataR = "TAATAT"
            tata2 = "TATAAA"
            tata2R = "AAATAT"
            threshold = 3
            for i in range(6):
                each_feature_vector.append(methods.string_matching(p[st:en], tata1, threshold))
                tata1 = tata1[1:] + tata1[:1]

                each_feature_vector.append(methods.string_matching(p[st:en], tataR, threshold))
                tataR = tataR[1:] + tataR[:1]

                each_feature_vector.append(methods.string_matching(p[st:en], tata2, threshold))
                tata2 = tata2[1:] + tata2[:1]

                each_feature_vector.append(methods.string_matching(p[st:en], tata2R, threshold))
                tata2R = tata2R[1:] + tata2R[:1]

            st = 20
            en = 35
            tata1 = "TTGACA"
            tataR = "ACAGTT"
            for i in range(6):
                each_feature_vector.append(methods.string_matching(p[st:en], tata1, threshold))
                tata1 = tata1[1:] + tata1[:1]

                each_feature_vector.append(methods.string_matching(p[st:en], tataR, threshold))
                tataR = tataR[1:] + tataR[:1]

            each_feature_vector.append(methods.string_matching(p[61:67], "AACGAT", threshold))


        # Feature 9
        # Position distance count of A,C,G,T
        # --------------------------------------------------------------------------------------
        if(feature_list[13]):
            each_feature_vector.append(methods.distance_count(p, 'A'))
            each_feature_vector.append(methods.distance_count(p[20:35], 'A'))
            each_feature_vector.append(methods.distance_count(p[45:55], 'A'))
            each_feature_vector.append(methods.distance_count(p, 'C'))
            each_feature_vector.append(methods.distance_count(p[20:35], 'C'))
            each_feature_vector.append(methods.distance_count(p[45:55], 'C'))
            each_feature_vector.append(methods.distance_count(p, 'G'))
            each_feature_vector.append(methods.distance_count(p[20:35], 'G'))
            each_feature_vector.append(methods.distance_count(p[45:55], 'G'))
            each_feature_vector.append(methods.distance_count(p, 'T'))
            each_feature_vector.append(methods.distance_count(p[20:35], 'T'))
            each_feature_vector.append(methods.distance_count(p[45:55], 'T'))


        # Feature 10
        # Dinucleotide Parameters Based on DNasel Digestion Data.
        # --------------------------------------------------------------------------------------
        if (feature_list[14]):
            each_feature_vector.append(methods.dinucleotide_value(p))
            each_feature_vector.append(methods.dinucleotide_value(p[20:35]))
            each_feature_vector.append(methods.dinucleotide_value(p[45:55]))


        newstr.append(each_feature_vector)

    return order, sequences, newstr