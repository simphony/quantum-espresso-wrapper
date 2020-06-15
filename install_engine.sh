./install_engine_requirements.sh

echo "Downloading a recent stable version"
wget https://github.com/QEF/q-e/releases/download/qe-6.5/qe-6.5-ReleasePack.tgz


echo "Extracting..."
tar  -xvzf qe-6.5-ReleasePack.tgz
cd qe-6.5
echo "Configuring..."
./configure
echo "Compiling..."
make all