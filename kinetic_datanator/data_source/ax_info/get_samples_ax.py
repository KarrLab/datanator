import requests
import os



def get_sample_text(ax_nums):
	for ax_num in ax_nums:
		print(ax_num)
		file = open("AX_samples_text/{}.txt".format(ax_num), "w")
		file.write(requests.get("https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/{}/samples".format(ax_num)).text)



if __name__ == '__main__':
	ax_nums = ["E-MTAB-2985",
		"E-MTAB-4980",
		"E-MTAB-5324",
		"E-MTAB-5181",
		"E-MTAB-4766",
		"E-MTAB-3226",
		"E-MTAB-4723",
		"E-MTAB-4722",]

	get_sample_text(ax_nums)
