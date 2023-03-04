#!/bin/bash

# dictionary of network-names & passwords accepted by each device....

declare -A networks=(
    ["network1_ssid"]="network1_password"
    ["network2_ssid"]="network2_password"
    ["network3_ssid"]="network3_password"
)

# Loop through all wireless network devices
for device in /sys/class/net/wlan*; do
    # Get device name
    device_name=$(basename $device)
    
    # Scan for available networks and get SSID names
    network_ssids=$(iw dev $device_name scan | grep "SSID:" | awk -F ":" '{print $2}')
    
    # The loop through available networks
    for ssid in $network_ssids; do
        # Check if network is in the list of networks to connect to
        if [ "${networks[$ssid]}" ]; then
            # Network SSID and password
            SSID=$ssid
            PASSWORD=${networks[$ssid]}
            # Configure wpa_supplicant (each device) with the corresponding config
            cat > /etc/wpa_supplicant/wpa_supplicant-$device_name.conf <<EOF
ctrl_interface=/run/wpa_supplicant
network={
        ssid="$SSID"
        psk="$PASSWORD"
}
EOF
            # Then connect to that network via said config
            wpa_supplicant -B -i $device_name -c /etc/wpa_supplicant/wpa_supplicant-$device_name.conf
            sleep 5
            # Then get an IP address from DHCP server
            dhclient $device_name
            # Once aquired, break the loop
            break
        fi
    done
done
