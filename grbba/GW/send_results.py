# -*- coding: UTF-8 -*-
'''
send alert email to BA
'''
import smtplib
from email.mime.text import MIMEText
mailto_list=["luoqi@ihep.ac.cn"]
mail_host="smtp.139.com" 
mail_user="luoqi_ba"    
mail_pass="luoqi@172"   
mail_postfix="139.com"

def send_mail(to_list,sub,content):
    me="<"+mail_user+"@"+mail_postfix+">"#"GW_uplim_REPORT"+
    msg = MIMEText(content,_subtype='plain',_charset='gb2312')
    msg['Subject'] = sub
    msg['From'] = me
    msg['To'] = ";".join(to_list)
    try:
        server = smtplib.SMTP()
        server.connect(mail_host)
        server.login(mail_user,mail_pass)
        server.sendmail(me, to_list, msg.as_string())
        server.close()
        return True
    except Exception, e:
        print str(e)
        return False
if __name__ == '__main__':
    if send_mail(mailto_list,"GW uplim Result","text to add"):
        print "send email successful"
    else:
        print "send email failed"
