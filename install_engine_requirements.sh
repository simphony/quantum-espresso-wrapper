echo "Installing engine requirements..."
platform = $(python3 -mplatform)

case $platform in
    *"Ubuntu"*)
        sudo apt-get install gfortran
        sudo apt-get install gcc
        sudo apt-get openmpi-bin
    ;;
#Support for more distributions to be added later
esac
