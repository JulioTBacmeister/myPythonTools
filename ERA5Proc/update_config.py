# Python Script: update_config.py
import yaml
from datetime import datetime, timedelta

def read_config_yaml(file_path):
    with open(file_path, "r") as config_file:
        config = yaml.safe_load(config_file)
    return config

def write_config_yaml(file_path, config):
    with open(file_path, "w") as config_file:
        yaml.dump(config, config_file)

def increment_month(config):
    # Increment the month and handle month/year change if needed
    # This is a simplistic implementation and does not handle all edge cases
    config['month'] += 1
    if config['month'] > 12:
        config['month'] = 1
        config['year'] += 1
    return config

def decrement_Resubmit(config):
    # Decrement the Resubmit counter
    config['Resubmit'] -= 1
    return config

def main():
    file_path = './config.txt'  # Specify the path to your config file
    config = read_config(file_path)
    config = increment_day(config)
    write_config(file_path, config)

if __name__ == "__main__":
    main()
