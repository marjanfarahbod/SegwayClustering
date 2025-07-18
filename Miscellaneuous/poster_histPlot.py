

# the hist data


baraxislabels = [5, 1, 15, 4, 6, 11, 12, 10, 2, 3, 14, 13, 0, 9, 7, 8]
barData = [.025, 0.12, .07, .06, .022, .06, .06, .002, .013, .11, .007, .06, .13, .06, .03, .17]


fix, axs = plt.subplots(1, figsize=(6,4))
#plt.grid()
plt.bar(range(len(baraxislabels)), barData)
#plt.ylim([15,110])
#plt.xticks(range(11))
axs.set_xticklabels(baraxislabels)
#axs.set_xticklabels(unique_track_list, rotation=90)
plt.tight_layout()
#plt.show()
figFile = '/Users/marjanfarahbod/Documents/talks/RSGDREAM2022/fig7.pdf'
print(figFile)
plt.savefig(figFile)
plt.close('all')
