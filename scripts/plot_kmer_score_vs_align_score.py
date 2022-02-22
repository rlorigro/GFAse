from matplotlib import pyplot
import numpy


def main():
    csv_path = "/home/ryan/data/test_gfase/paolo_ul_run9/phased_k31_double_coverage_renamed/unphased_parental_counts_aligned.csv"

    fig = pyplot.figure()
    axes = pyplot.axes()

    x = list()
    y_score = list()
    y_unique_score = list()
    sizes = list()
    names = list()

    # 0 name
    # 1 maternal_count
    # 2 paternal_count
    # 3 unique_maternal_count
    # 4 unique_paternal_count
    # 5 score
    # 6 unique_score
    # 7 color
    # 8 ref_name
    # 9 ref_score
    # 10 length
    with open(csv_path, 'r') as file:
        for l,line in enumerate(file):
            if l == 0:
                continue

            data = line.strip().split(',')

            ref_name = data[8]

            # if "X" not in ref_name and "Y" not in ref_name:
            #     continue
            # if "X" not in ref_name:
            #     continue
            # if "Y" not in ref_name:
            #     continue

            maternal_count = int(float(data[1]))
            paternal_count = int(float(data[2]))

            if maternal_count + paternal_count < 30:
                continue

            print(line.strip())

            kmer_score = float(data[5])
            unique_kmer_score = float(data[6])
            alignment_score = float(data[9])
            query_length = float(data[10])

            print(maternal_count, paternal_count, kmer_score, unique_kmer_score, alignment_score, query_length)

            x.append(alignment_score)
            y_score.append(kmer_score)
            y_unique_score.append(unique_kmer_score)
            sizes.append(query_length)
            names.append(ref_name)

            print()

    # y = y_score
    y = y_unique_score

    threshold = 0.02

    log_sizes = [max(1,numpy.log10(x)-2)**3 for x in sizes]

    pyplot.scatter(x, y, marker='o', s=log_sizes)

    axes.axhline(0, linestyle='--', linewidth=0.6, color='gray')
    axes.axvline(0, linestyle='--', linewidth=0.6, color='gray')

    total_unphased = 0
    total_paternal = 0
    total_maternal = 0

    for i in range(len(names)):
        name = names[i]

        if "X" in name:
            # if sizes[i] > 50000:
            #     pyplot.text(x[i], y[i], "X %.1f" % (float(sizes[i])/1000))
            pass

        elif "Y" in name:
            if sizes[i] > 50000:
                pyplot.text(x[i], y[i], "Y %.1f" % (float(sizes[i])/1000))
            pass

        # elif sizes[i] > 100000:
        #     pyplot.text(x[i], y[i], "%s %.1f" % (name, float(sizes[i])/1000))

        if sizes[i] > 100_000:
            if y[i] > threshold:
                total_paternal += sizes[i]
            elif y[i] < -threshold:
                total_maternal += sizes[i]
            else:
                total_unphased += sizes[i]
        else:
            total_unphased += sizes[i]

    print("total_unphased", total_unphased)
    print("total_paternal", total_paternal)
    print("total_maternal", total_maternal)

    fig2 = pyplot.figure()
    axes2 = pyplot.axes()

    total = sum(sizes)

    x_prev = 0
    s_prev = None
    for s in sorted(sizes, reverse=True):
        x0 = x_prev
        x1 = x_prev + s
        pyplot.plot([x0,x1],[s,s], color="C1")

        if s_prev is not None:
            pyplot.plot([x0,x0],[s_prev,s], color="C1")

        x_prev = x1
        s_prev = s

    axes2.axvline(float(total)/2.0, linestyle='--', linewidth=0.6, color='gray')
    pyplot.show()
    pyplot.close()


if __name__ == "__main__":
    main()
