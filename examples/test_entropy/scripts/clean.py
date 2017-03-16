#
#
#

# import my stuff
from common_batch import delete_file
from common_batch import delete_directory_content

# delete created folders and files
delete_directory_content("batch/")
delete_directory_content("batch_out/")
delete_directory_content("results/")
delete_directory_content("log/")
delete_directory_content("shortinfo/")
delete_file("shortinfo_final.txt")



