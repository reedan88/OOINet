# OOINet
The modules and tools within this repo are designed to assist in requesting, importing, downloading, and visualizing data from the Ocean Observatories Initiative API by M2M requests.

The Ocean Observatories Initiative (OOI) is a science-driven ocean observing network that delivers real-time data from more than 800 instruments to address critical science questions regarding the world’s ocean. OOI data are freely available online to anyone with an Internet connection. OOI features two methods for accessing data: via the web portal at [ooinet.oceanobservatories.org ](https://ooinet.oceanobservatories.org) or via OOINet's [API](https://oceanobservatories.org/ooi-m2m-interface/).

<figure>
<img src="figures/overview.png">
</figure>

The ooinet tool is designed to simplify and streamline the searching, requesting, and downloading of netCDF datasets from OOINet via the API. It is designed for interactive computing in jupyter notebooks for the science end-user. It allows for searching of specific datasets based on the location and/or specific instrument and built-in functions reprocess datasets to allow for multiple overlapping netCDF files to be loaded remotely into an xarray dataset. It will also retrieve associated metadata and vocabulary for specific datasets.

---
## Setup
1. **OOI Credentials**

Before using the tool you will need to register an account on [ooinet.oceanobservatories.org ](https://ooinet.oceanobservatories.org) and get an API username and token that allow. Then, under your **User Profile** you will be assigned an **API Username** and **API Token.** If you have not done so, follow these steps:
  * Navigate to the [OOI Data Portal](ooinet.oceanobseravtories.org)
  * At the upper right of the page, click the "login"
  * Once you are logged in, navigate to your **User Profile** in the upper right of the screen

  <figure>
  <img src="figures/user_profile_screen.png">
  </figure><b>

  * On your user profile page, save the "API Username" and "API Token"

  <figure>
  <img src="figures/api_credentials.png">
  </figure>

2. **Save your credetials**

  Following the instructions from the ooi-data-explorations home page, we'll utilize the ```netrc``` module to store your credentials on your local machine.

  ```
  cd ~
  touch .netrc
  chmod 600 .netrc
  cat <<EOT >> .netrc
  machine ooinet.oceanobservatories.org
      login <API Username>
      password <API Token>
  EOT
  ```
---
## Usage

An [example jupyter notebook](examples/) is provided in the examples subfolder.

The key parameters which the OOINet API requires is the "reference designator." A reference designator may be thought of as a type of instrument located at a fixed location and depth. In the example notebook, we use the CP01CNSM-RID27-03-CTDBPC000, which is the CTD located at 7 meters depth on the Pioneer Array Central Surface Mooring at approximately (latitude, longitude) of (40.14, -70.7783).


## Meta
Andrew Reed <br>
mailto:areed@whoi.edu <br>
Distributed under the ```GPU GPLv3``` license. See ```setup.py``` for more information. <br>
[https://github.com/reedan88/OOINet](https://github.com/reedan88/OOINet)


## Contributing
Contributions are always welcome, particularly assistance with packaging for distribution via PyPI.
1. Fork it
2. Create your feature branch (```git checkout -b feature/```)
3. Commit your changes (```git commit -m "Added a feature"```)
4. Push the branch to your fork (```git push origin feature/```)
5. Create a pull request
