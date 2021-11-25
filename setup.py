from email.utils import getaddresses

import toml
from setuptools import find_packages, setup


module = toml.load("pyproject.toml")
module_tool_poetry = module["tool"]["poetry"]
module_name = module_tool_poetry["name"]

module_package_from = module_tool_poetry["packages"][0]["from"]

module_package_include = module_tool_poetry["include"]

module_install_requires = [
    k + v.replace("^", "!=").replace("*", "")
    if type(v) == str
    else k + v["version"].replace("^", "!=").replace("*", "")
    for (k, v) in module_tool_poetry["dependencies"].items()
    if k != "python"
]

module_devel_requires = [
    k + v.replace("^", "!=").replace("*", "")
    if type(v) == str
    else k + v["version"].replace("^", "!=").replace("*", "")
    for (k, v) in module_tool_poetry["dev-dependencies"].items()
]


setup(
    name=module_name,
    version=module_tool_poetry["version"],
    author=", ".join(e[0] for e in getaddresses(module_tool_poetry["authors"])),
    author_email=", ".join(e[1] for e in getaddresses(module_tool_poetry["authors"])),
    license=module_tool_poetry["license"],
    description=module_tool_poetry["description"],
    long_description=open(module_tool_poetry["readme"]).read(),
    platforms="all",
    classifiers=module_tool_poetry["classifiers"],
    package_dir={"": module_package_from},
    packages=[module_tool_poetry["packages"][0]["include"]],
    include_package_data=True,
    package_data={"": module_package_include},
    keywords=", ".join(module_tool_poetry["keywords"]),
    python_requires=">=3.7, <4",
    install_requires=module_install_requires,
    test_suite="tests",
    extras_require={
        "develop": module_devel_requires,
    },
    project_urls={
        "Documentation": module_tool_poetry["homepage"],
        "Source": module_tool_poetry["repository"],
    },
)
