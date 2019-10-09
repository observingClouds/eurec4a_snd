"""
Script to create or adjust necessary config file
"""

import configparser


def read_config(cfg_file):
    """
    read config file
    """
    config = configparser.ConfigParser()
    config.read(cfg_file)
    return config


def rqst_user_input(property, current_value):
    """
    Ask user to modify the property
    given and retuns the modification
    """
    print("Currently: {prop} = {value}".
          format(prop=property,
                 value=current_value))
    update = input("Please modify (leave blank to leave unchanged):")
    if update == '':
        return current_value
    else:
        return update


def adjust_config(cfg_instance):
    """
    Modification of configurations
    by user input
    """
    cfg_sections = cfg_instance.sections()
    for section in cfg_sections:
        for field in cfg_instance[section]:
            update = rqst_user_input(field, cfg_instance[section][field])
            cfg_instance[section][field] = update
    return cfg_instance


def write_config(cfg_instance, cfg_file):
    """
    Write updated configuration back to config
    file
    """
    print("Writing updated config file to {}".format(cfg_file))
    with open(cfg_file, 'w') as configfile:
        cfg_instance.write(configfile)


def update_config(cfg_file, cfg_file_new_location):
    """
    Update config file or create a personal copy
    in the user space
    """
    current_config = read_config(cfg_file)
    updated_config = adjust_config(current_config)
    write_config(updated_config, cfg_file_new_location)


if __name__ == '__main__':
    pass
