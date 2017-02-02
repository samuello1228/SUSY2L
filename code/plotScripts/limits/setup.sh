if [ ! -e config/HistFactorySchema.dtd ]; then
	echo 'Linking HistFactorySchema.dtd from HistFitter'
	if [ ! -d config ]; then mkdir config; fi
	ln -s ../../../HistFitter/config/HistFactorySchema.dtd config
else
	echo "HistFactorySchema.dtd already exist. Good!"
fi
