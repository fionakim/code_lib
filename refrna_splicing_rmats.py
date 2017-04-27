# -*- coding: utf-8 -*-
# __author__ = fiona
# time: 2017/4/27 16:31

import re, os, Bio, argparse, sys, fileinput, urllib2
from biocluster.api.database.base import Base, report_check
import os
import datetime
import types
from biocluster.config import Config
from bson.son import SON
from bson.objectid import ObjectId
import bson.binary
from cStringIO import StringIO
import re
import json
import pandas as pd
import numpy as np
from collections import Counter


class RefrnaSplicing(Base):
    def __init__(self, bind_object):
        super(RefrnaSplicing, self).__init__(bind_object)
        self._db_name = Config().MONGODB + '_ref_rna'
    
    @report_check
    def add_splicing(self, params=None, rmats_dir=None, name=None, samples=None, major=True):
        task_id = self.bind_object.sheet.id
        project_sn = self.bind_object.sheet.project_sn
        insert_data = {
            'project_sn': project_sn,
            'task_id': task_id,
            'name': name if name else 'express_matrix_' + str(datetime.datetime.now().strftime("%Y%m%d_%H%M%S")),
            'desc': '可变剪接结果主表',
            'created_ts': datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'params': (
                json.dumps(params, sort_keys=True, separators=(',', ':')) if isinstance(params, dict) else params),
            'specimen': samples,
            'status': 'end',
        }
        collection = self.db['sg_splicing_rmats']
        splicing_rmats_id = collection.insert_one(insert_data).inserted_id
        if major:
            rmats_event_files = [f for f in os.listdir(os.path.join(rmats_dir, 'ASEvents')) if
                                 re.match(r'fromGTF\.(novelEvents)?\.(A5SS|A3SS|SE|RI|MXE)\.alter_id\.txt', f)]
            rmats_mats_files = [f for f in os.listdir(os.path.join(rmats_dir, 'MATS_output')) if
                                re.match(
                                    r'(A5SS|A3SS|SE|RI|MXE)\.MATS\.(ReadsOnTargetAndJunctionCounts|JunctionCountOnly)\.alter_id\.psi_info\.txt',
                                    f)]
            
        
        return splicing_rmats_id
    
    def add_splicing_detail(self, splicing_id):
        if not isinstance(splicing_id, ObjectId):
            if isinstance(splicing_id, types.StringTypes):
                splicing_id = ObjectId(splicing_id)
            else:
                raise Exception('splicing_id必须为ObjectId对象或其对应的字符串！')
    
    def add_splicing_stats(self, splicing_id):
        if not isinstance(splicing_id, ObjectId):
            if isinstance(splicing_id, types.StringTypes):
                splicing_id = ObjectId(splicing_id)
            else:
                raise Exception('splicing_id必须为ObjectId对象或其对应的字符串！')
    
    def add_splicing_psi(self, splicing_id):
        if not isinstance(splicing_id, ObjectId):
            if isinstance(splicing_id, types.StringTypes):
                splicing_id = ObjectId(splicing_id)
            else:
                raise Exception('splicing_id必须为ObjectId对象或其对应的字符串！')
