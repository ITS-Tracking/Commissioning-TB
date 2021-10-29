# make sure to run this in a directory named after the run you want to analyse

mkdir $1
ln -s ~/alice/run3/data/o2sim_geometry.root $1/
ln -s ~/alice/run3/data/o2sim_grp2.root $1/o2sim_grp.root
ln -s ~/alice/run3/data/ctf_dictionary.root $1/
ln -s ~/alice/run3/data/ITSdictionary.bin $1/
ln -s ~/alice/run3/data/matbud.root $1/
ln -s $(pwd)/command.sh $1
xrdfs root://eosaliceo2 ls  /eos/aliceo2/raw/2021/OCT/${PWD##*/}/raw/$1 | grep ctf > $1/data_all.lst
cd $1
python3 ../cleanup.py > data.lst
cd ..
