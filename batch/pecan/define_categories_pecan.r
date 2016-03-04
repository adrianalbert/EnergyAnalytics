# define_categories_pecan.r
#
# Defines appliance categories for the Pecan St dataset.
#
# Adrian Albert
# Last modified: May 2014.
# ---------------------------------------------------------

# __________________________________________________
# Define appliance categories

# some interesting components
# appliances   = as.character(read.csv('~/energy-data/pecan_street/metadata/appliances.csv')$Appliance)
select.keep  = c('dataid', 'localminute', 'use')
select.AC    = c("air1", "air2", "air3", "housefan1")
# select.HV    = c("furnace1", "furnace2", "heater1")
select.HV    = c("heater1", "airwindowunit1", "furnace1", "furnace2")
select.light = c("lights_plugs1", "lights_plugs2", "lights_plugs3", "lights_plugs4", "lights_plugs5", "lights_plugs6",
                 "outsidelights_plugs1", "outsidelights_plugs2")
select.alwOn = c('refridgerator1', 'refridgerator2', 'winecooler1', 'aquarium1',
                 "freezer1")
select.sched = c("pool1", "pool2", 'sprinkler1', "poolpump1", "pump1")
select.total = c('use')
select.dhw   = c('waterheater1', 'waterheater2')
select.user  = c("bathroom1", "bathroom2", "bedroom1", "bedroom2", "bedroom3", "bedroom4", "bedroom5",
                 "clotheswasher1", "clotheswasher_dryg1", "diningroom1", "diningroom2", "dishwasher1",
                 "disposal1", "drye1", "dryg1", "garage1", "garage2", "icemaker1", "jacuzzi1", 
                 "kitchenapp1", "kitchenapp2", "livingroom1", "livingroom2", "heater1",
                 "microwave1", "office1", "oven1", "poollight1",  "range1", "security1", "shed1", "utilityroom1", "venthood1")
select.solar = c('gen')
select.ev    = c('car1')
