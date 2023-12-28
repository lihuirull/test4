# -*- coding:UTF-8 -*-

# def refreshStandardName(predictedProteinType,dic): #更新了字典dic和列表predictedProteinType
#     predictedProteinTypeNew = []
#     for eachType in predictedProteinType:
#         #更新了字典dic和predictedProteinType，改变了NS翻译出两种蛋白NS1和NS2，MP对应的M1M2的名字，
#         # 例如MP的querySeq1变为querySeq_M1和querySeq_M2两个，蛋白质类型也由MP改为M1和M2
#         if eachType[1] == "NS":
#             predictedProteinTypeNew.append((eachType[0] + "_NS1", "NS1"))
#             predictedProteinTypeNew.append((eachType[0] + "_NS2", "NS2"))
#             dic[">" + eachType[0] + "_NS1"] = dic[">" + eachType[0]]
#             dic[">" + eachType[0] + "_NS2"] = dic[">" + eachType[0]]
#             dic.pop(">" + eachType[0])
#
#         elif eachType[1] == "MP":
#             predictedProteinTypeNew.append((eachType[0] + "_M1", "M1"))
#             predictedProteinTypeNew.append((eachType[0] + "_M2", "M2"))
#             dic[">" + eachType[0] + "_M1"] = dic[">" + eachType[0]]
#             dic[">" + eachType[0] + "_M2"] = dic[">" + eachType[0]]
#             dic.pop(">" + eachType[0])
#         elif eachType[1] == "PB1":
#             predictedProteinTypeNew.append((eachType[0] + "_PB1", "PB1"))
#             predictedProteinTypeNew.append((eachType[0] + "_PB1-F2", "PB1-F2"))
#             dic[">" + eachType[0] + "_PB1"] = dic[">" + eachType[0]]
#             dic[">" + eachType[0] + "_PB1-F2"] = dic[">" + eachType[0]]
#             dic.pop(">" + eachType[0])
#         elif eachType[1] == "PA":
#             predictedProteinTypeNew.append((eachType[0] + "_PA", "PA"))
#             predictedProteinTypeNew.append((eachType[0] + "_PA-X", "PA-X"))
#             dic[">" + eachType[0] + "_PA"] = dic[">" + eachType[0]]
#             dic[">" + eachType[0] + "_PA-X"] = dic[">" + eachType[0]]
#             dic.pop(">" + eachType[0])
#         else:
#             predictedProteinTypeNew.append(eachType)
#     return predictedProteinTypeNew #dic因为是引用，所以无需return出函数，已经在函数发生了改变

def refreshStandardName(predictedProteinType, dic):
    """
    根据预测的蛋白类型更新蛋白质名称列表和字典。

    对于特定的蛋白类型（NS, MP, PB1, PA），函数会将每个蛋白质名称分割为更具体的子类型（如NS1, NS2等）。
    同时，它也会在字典中更新这些蛋白质的记录。

    :param predictedProteinType: 被预测的蛋白类型的列表，列表中的每个元素是一个包含蛋白质名称和类型的元组。
    :param dic: 一个字典，键是蛋白质的名称，值是相关的数据。
    :return: 更新后的蛋白质类型列表。
    """

    predictedProteinTypeNew = []
    for eachType in predictedProteinType:
        # 根据蛋白质类型进行不同的处理
        if eachType[1] == "NS":
            # NS类型的蛋白质分成NS1和NS2
            predictedProteinTypeNew.append((eachType[0] + "_NS1", "NS1"))
            predictedProteinTypeNew.append((eachType[0] + "_NS2", "NS2"))
            dic[">" + eachType[0] + "_NS1"] = dic[">" + eachType[0]]
            dic[">" + eachType[0] + "_NS2"] = dic[">" + eachType[0]]
            dic.pop(">" + eachType[0])

        elif eachType[1] == "MP":
            # MP类型的蛋白质分成M1和M2
            predictedProteinTypeNew.append((eachType[0] + "_M1", "M1"))
            predictedProteinTypeNew.append((eachType[0] + "_M2", "M2"))
            dic[">" + eachType[0] + "_M1"] = dic[">" + eachType[0]]
            dic[">" + eachType[0] + "_M2"] = dic[">" + eachType[0]]
            dic.pop(">" + eachType[0])

        elif eachType[1] == "PB1":
            # PB1类型的蛋白质分成PB1和PB1-F2
            predictedProteinTypeNew.append((eachType[0] + "_PB1", "PB1"))
            predictedProteinTypeNew.append((eachType[0] + "_PB1-F2", "PB1-F2"))
            dic[">" + eachType[0] + "_PB1"] = dic[">" + eachType[0]]
            dic[">" + eachType[0] + "_PB1-F2"] = dic[">" + eachType[0]]
            dic.pop(">" + eachType[0])

        elif eachType[1] == "PA":
            # PA类型的蛋白质分成PA和PA-X
            predictedProteinTypeNew.append((eachType[0] + "_PA", "PA"))
            predictedProteinTypeNew.append((eachType[0] + "_PA-X", "PA-X"))
            dic[">" + eachType[0] + "_PA"] = dic[">" + eachType[0]]
            dic[">" + eachType[0] + "_PA-X"] = dic[">" + eachType[0]]
            dic.pop(">" + eachType[0])

        else:
            # 其他类型不变更
            predictedProteinTypeNew.append(eachType)

    # 返回更新后的蛋白质类型列表，dic作为可变对象已在函数中直接修改。
    return predictedProteinTypeNew
