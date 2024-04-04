import matplotlib.pyplot as plt


def plot_promoter_overlap(df_promoters, streptococcus_promoter, extension_bp, save_directory, index):
    # 确保保存目录存在
    if not os.path.exists(save_directory):
        os.makedirs(save_directory)

    # 获取单个预测启动子
    predicted = df_promoters.loc[index]
    pred_start, pred_end = predicted['start'], predicted['end']

    plt.figure(figsize=(10, 2))
    plt.title(f"Predicted vs Real Promoters (Index: {index})")
    found_overlap = False

    # 遍历所有真实启动子，查找重叠
    for real_index, real in streptococcus_promoter.iterrows():
        real_start, real_end = real['Promoter Start'] - extension_bp, real['Promoter End'] + extension_bp

        if pred_end >= real_start and pred_start <= real_end:
            # 绘制重叠的启动子
            plt.plot([pred_start, pred_end], [1, 1], 'r-', linewidth=5, label='Predicted Promoter')
            plt.plot([real_start, real_end], [1.1, 1.1], 'g-', linewidth=5, label='Real Promoter')
            found_overlap = True
            break  # 假设每个预测启动子只与一个真实启动子匹配

    if found_overlap:
        plt.legend()
        plt.ylim(0.9, 1.2)
        plt.xlim(min(pred_start, real_start) - 10, max(pred_end, real_end) + 10)
        plt.yticks([])
        plt.xlabel('Genome Position')
        plt.show()
        # 构建文件名和保存路径
        file_name = f"promoter_overlap_{index}.png"
        save_path = os.path.join(save_directory, file_name)

        # 保存并关闭图表
        plt.savefig(save_path)
        plt.close()
    else:
        plt.close()  # 如果没有找到重叠，关闭图表而不保存


def analyze_predictions(df_promoters, streptococcus_promoter, extension_bp):
    used_real_indices = set()
    matched_predictions = 0
    total_predicted_bases = 0
    true_positives = 0
    false_positives = 0

    for index, predicted in df_promoters.iterrows():
        pred_start, pred_end = predicted['start'], predicted['end']
        total_predicted_bases += pred_end - pred_start + 1

        for real_index, real in streptococcus_promoter.iterrows():
            if real_index in used_real_indices:
                continue

            real_start, real_end = real['Promoter Start'] - extension_bp, real['Promoter End'] + extension_bp

            overlap_start = max(pred_start, real_start)
            overlap_end = min(pred_end, real_end)
            if overlap_end >= overlap_start:
                matched_predictions += 1
                used_real_indices.add(real_index)
                
                overlap_length = overlap_end - overlap_start + 1
                true_positives += overlap_length
                break  # Once a match is found, we break the loop

    false_positives = total_predicted_bases - true_positives
    return matched_predictions, total_predicted_bases, true_positives, false_positives
def run_analysis(df_promoters, streptococcus_promoter, extension_bp, save_directory):
    # Analyze predictions
    matched, total, tp, fp = analyze_predictions(df_promoters, streptococcus_promoter, extension_bp)
    
    # Iterate over each promoter and plot the overlaps
    for index in df_promoters.index:
        plot_promoter_overlap(df_promoters, streptococcus_promoter, extension_bp, save_directory, index)

    # Print analysis results
    print(f"Matched Predictions: {matched}")
    print(f"Total Predicted Bases: {total}")
    print(f"True Positives: {tp}")
    print(f"False Positives: {fp}")
#streptococcus_promoter = pd.read_csv('/home/hanzequan/test_bectiral/bxb1_tfbs.csv')
#streptococcus_promoter.rename(columns={
#    'sequence': 'Promoter Sequence',
#    'start': 'Promoter Start',
#    'end': 'Promoter End'
#}, inplace=True)


#run_analysis(df_promoters, meme_results, extension_bp=0, save_directory='/home/hanzequan/test_bectiral/EJ_1')
