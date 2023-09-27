#!/usr/bin/env python
# coding: utf-8
import sys
import os
from utils import *


def main(config):
    incoming_data_path = config['paths']['incoming_data']
    processed_data_path = config['paths']['data']
    


if __name__ == '__main__':
    CONFIG = load_config()
    main(CONFIG)