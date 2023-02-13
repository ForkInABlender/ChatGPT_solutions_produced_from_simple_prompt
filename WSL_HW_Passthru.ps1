$wslConfigFile = "C:\etc\wsl.conf"

# Get a list of all devices connected to the system
$deviceList = Get-WmiObject -Class Win32_PnPEntity | Select-Object -Property DeviceID, PNPDeviceID

# Write the WSL2 configuration file
$config = "[wsl2]"
$config | Out-File -FilePath $wslConfigFile

# For each device in the list, generate a configuration line
foreach ($device in $deviceList) {
  $configLine = "device=$($device.DeviceID):/dev/$($device.PNPDeviceID)"
  $configLine | Out-File -FilePath $wslConfigFile -Append
}
