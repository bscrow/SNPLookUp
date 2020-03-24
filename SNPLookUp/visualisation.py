import matplotlib.pyplot as plt
import os


def visualise(query_region, save_to_jpeg=False, filename="visualisation"):
    if type(query_region) == dict:
        temp = []
        for qreg in query_region.values():
            temp.extend(qreg)
        query_region = temp
    height = [len(qreg.get_snps()) for qreg in query_region]
    labels = [qreg.label for qreg in query_region]
    xpos = [a * 2 + 1 for a in range(len(height))]
    bar_width = 0.9
    plt.figure(figsize=(min(len(query_region) * 2, 10), 10))
    plt.style.use('seaborn-pastel')

    plt.bar(xpos, height, width=bar_width)

    plt.title(f"Visualisation of query region\n", fontsize=20)

    for i in range(len(height)):
        plt.text(x=xpos[i], y=height[i] + 2, s=str(height[i]), size=16, horizontalalignment='center')
    plt.xticks(xpos, labels, fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel("Query Region Label", fontsize=16)
    plt.ylabel("Number of SNPs", fontsize=16)
    plt.grid()

    if save_to_jpeg:
        if not os.path.exists(os.path.join("..", "output", "visualisation")):
            os.mkdir(os.path.join("..", "output", "visualisation"))
        plt.savefig(os.path.join("..", "output", "visualisation", filename), dpi=400)

    plt.show()
