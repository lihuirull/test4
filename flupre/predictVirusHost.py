def getVirusHost(predictedProteinType,Host):
    from collections import Counter
    HostNew = []
    for eachProtein in predictedProteinType:
        # if eachProtein[1] != "M2" and eachProtein[1] != "PA-X" and eachProtein[1] != "PB1-F2" and eachProtein[1] != "NS1" and eachProtein[1] != "NS2"and eachProtein[1] != "M1"and eachProtein[1].strip() != ""and eachProtein[1] != "Unknown":
        if eachProtein[1] not in ['M2',"PA-X","PB1-F2","NS1","NS2","M1","Unknown",'']:
        #去掉短蛋白，防止防止短蛋白造成结论错误
            for eachHost in Host:
                if eachHost[0] == eachProtein[0]:
                    HostNew.append(eachHost)
    if len(HostNew) < 1:
        HostNew = Host
    virusHost = Counter(host[1] for host in HostNew).most_common()[:]  # 该病毒中所有蛋白预测的宿主，出现频路最高的那个选为病毒的宿主
    sum = 0
    for eachHostFrequency in virusHost:
        sum = float(eachHostFrequency[1])+sum
    virusHostNew = []
    n = 0
    for eachHostFrequency in virusHost:
        HostType = eachHostFrequency[0]
        probability = '%.3f' % (eachHostFrequency[1]/sum)
        virusHostNew.append((HostType,str(probability)))
        n = n+1
        if n>=3:
            break
    return str(virusHostNew).strip('[]'),HostNew
    
    
def getVirusHostLowLevel(predictedProteinType,Host):
    from collections import Counter
    HostNew = []
    for eachProtein in predictedProteinType:
        # if eachProtein[1] != "M2" and eachProtein[1] != "PA-X" and eachProtein[1] != "PB1-F2" and eachProtein[1] != "NS1" and eachProtein[1] != "NS2"and eachProtein[1] != "M1"and eachProtein[1].strip() != ""and eachProtein[1] != "Unknown":
        if eachProtein[1] not in ['M2',"PA-X","PB1-F2","NS1","NS2","M1","Unknown",'']:
        #去掉短蛋白，防止防止短蛋白造成结论错误
            for eachHost in Host:
                if eachHost[0] == eachProtein[0]:
                    HostNew.append(eachHost)
    if len(HostNew) < 1:
        HostNew = Host
    virusHost = Counter(host[1] for host in HostNew).most_common()[:]  # 该病毒中所有蛋白预测的宿主，出现频路最高的那个选为病毒的宿主
    sum = 0
    for eachHostFrequency in virusHost:
        sum = float(eachHostFrequency[1])+sum
    virusHostNew = []
    n = 0
    for eachHostFrequency in virusHost:
        HostType = eachHostFrequency[0]
        probability = '%.3f' % (eachHostFrequency[1]/sum)
        virusHostNew.append((HostType,str(probability)))
        n = n+1
        if n>=30:
            break
    return str(virusHostNew).strip('[]'),HostNew
