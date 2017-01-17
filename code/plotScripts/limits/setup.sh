if [ ! config/HistFactorySchema.dtd ]; then
	echo 'Linking HistFactorySchema.dtd from HistFitter'
# 	ln -s ../../../HistFitter/config/HistFactorySchema.dtd .
else
	echo "HistFactorySchema.dtd already exist. Good!"
fi
