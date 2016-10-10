#
#
#

# import my stuff
from myscripts import delete_file
from myscripts import delete_directory_content

# delete created folders and files
print "- delete folder batch/"
delete_directory_content("batch/")

print "- delete folder batch_out/"
delete_directory_content("batch_out/")

print "- delete folder results/"
delete_directory_content("results/")

print "- delete folder log/"
delete_directory_content("log/")

print "- delete folder shortinfo/"
delete_directory_content("shortinfo/")

print "- delete shortinfo.txt"
delete_file("shortinfo_final.txt")



