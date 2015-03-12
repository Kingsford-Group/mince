echo -e "Preparing binary release\n=====================\n"

# create the binary directory
mkdir $HOME/Mince-latest_ubuntu-12.04
mkdir $HOME/Mince-latest_ubuntu-12.04/bin
mkdir $HOME/Mince-latest_ubuntu-12.04/lib

# copy over the executable 
echo -e "Copying over the binary and Intel TBB libraries\n"
cp $TRAVIS_BUILD_DIR/build/srcmince $HOME/Mince-latest_ubuntu-12.04/bin/
cp $TRAVIS_BUILD_DIR/lib/libtbb* $HOME/Mince-latest_ubuntu-12.04/lib/

# copy other dependencies (shared libraries)
echo -e "Copying over other shared library dependencies\n"
bash $TRAVIS_BUILD_DIR/scripts/cpld.bash $TRAVIS_BUILD_DIR/build/src/mince $HOME/Mince-latest_ubuntu-12.04/lib/
echo -e "Removing dangerous dependencies\n"
rm $HOME/Mince-latest_ubuntu-12.04/lib/libc.so.6
rm $HOME/Mince-latest_ubuntu-12.04/lib/ld-linux-x86-64.so.2
rm $HOME/Mince-latest_ubuntu-12.04/lib/libdl.so.2

# now make the tarball
echo -e "Making the tarball\n"
cd $HOME
tar czvf Mince-latest_ubuntu-12.04.tar.gz Mince-latest_ubuntu-12.04

echo -e "Pushing the tarball to GitHub\n"
# Since it's currently unclear to me how to overwrite an asset via the GitHub
# API, the following code deletes the old asset, and uploads the new one in its place
# Get the previous asset id of the tarball
echo -e "Getting previous asset ID\n"
ASSETID=`curl -s -X GET https://api.github.com/repos/Kingsford-Group/mince/releases/1043508/assets | grep "\"id" | head -1 | awk '{gsub(/,$/,""); print $2}'`
# Delete the previous tarball
echo -e "Deleting previous asset\n"
curl -X DELETE -H "Authorization: token ${MINCE_PUSH_KEY}" https://api.github.com/repos/Kingsford-Group/mince/releases/assets/$ASSETID
# Upload the new tarball
echo -e "Uploading new asset\n"
curl -X POST --data-binary "@Mince-latest_ubuntu-12.04.tar.gz" https://uploads.github.com/repos/Kingsford-Groupmince/releases/1043508/assets?name=Mince-latest_ubuntu-12.04.tar.gz --header "Content-Type:application/gzip" -H "Authorization: token ${MINCE_PUSH_KEY}"
echo -e "Done!\n"

