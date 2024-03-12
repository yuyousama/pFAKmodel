# -*- coding: utf-8 -*-
"""
@Time    : 2022/10/29 17:05
@Author  : Qingdoors
"""
# 新模型
import numpy as np
import pylab as pl
import time


def stable(data, delta_t, time=300, easy_mode=False):
    threhold = 0.001
    if easy_mode:
        time /= 10
        threhold *= 50
    num_of_point = round(time / delta_t)
    if len(data) > num_of_point:
        if (max(data[-num_of_point:]) - min(data[-num_of_point:])) / min(data[-num_of_point:]) <= threhold:
            return True
    return False

class model:
    def __init__(self, delta_t, sets_of_k, show=True, pre_set=0):
        '''
        :param delta_t: 时间间隔
        :param sets_of_k:k值，单位 s-1？
        k0 off_pFAK的kin
        k1 on_pFAK的kout
        k2 off_FAK的kin，基础值
        k3 AP影响off_FAK的kin
        k4 on_FAK的kout
        k5 off_pFAK的去磷酸化
        k6 on_FAK的磷酸化，基础值
        k7 AP影响on_FAK的磷酸化
        k8 崩解率，0到0.2（0.25？）
        k9 k2的变化幅度
        k10 k2的变化速度
        :param show:展示细节
        :param preset:预设FAK量，减少计算量
        '''
        self.delta_t = delta_t
        self.sets_of_k = sets_of_k
        self.pre_set = pre_set
        self.show = show
        self.statu = 'dark'
        self.target_k2 = self.sets_of_k[2]
        self.p = self.sets_of_k[8]
        if show:
            print('Model initializing...')
        self.get_stable()

        if show:
            print('Model Ready')

    def get_stable(self):
        if np.sum(self.pre_set) < 500:
            self.off_pFAK = 250
            self.on_pFAK = 250
            self.on_FAK = 250
            self.off_FAK = 250
        else:
            if self.show:
                print('Loading...')
            self.off_pFAK = self.pre_set[0]
            self.on_pFAK = self.pre_set[1]
            self.on_FAK = self.pre_set[2]
            self.off_FAK = self.pre_set[3]
        time = 0

        t = []
        off_pFAK_list = []
        on_pFAK_list = []
        on_FAK_list = []
        off_FAK_list = []

        while True:
            if self.show:
                print('\r', time, end='', flush=True)
            off_pFAK_list.append(self.off_pFAK)
            on_pFAK_list.append(self.on_pFAK)
            on_FAK_list.append(self.on_FAK)
            off_FAK_list.append(self.off_FAK)
            t.append(time)
            self.update()
            time += self.delta_t
            if stable(off_pFAK_list, self.delta_t) and stable(on_pFAK_list, self.delta_t) \
                    and stable(on_FAK_list,self.delta_t) and stable(off_FAK_list, self.delta_t):
                break
            if min(self.off_pFAK, self.on_pFAK, self.on_FAK, self.off_FAK) < 0:
                if self.show:
                    print()
                    print('Too Low!!!')
                    print(self.off_pFAK, self.on_pFAK, self.on_FAK, self.off_FAK)
                    pl.plot(t, off_pFAK_list, label='off_pFAK', linewidth=3, c='#1f78b4')
                    pl.plot(t, on_pFAK_list, label='on_pFAK', linewidth=3, c='#33a02c')
                    pl.plot(t, on_FAK_list, label='on_FAK', linewidth=3, c='#e31a1c')
                    pl.plot(t, off_FAK_list, label='off_FAK', linewidth=3, c='#6a3d9a')
                    pl.legend(loc='best')
                    pl.show()

        if self.show:
            print()
            print('平衡态FAK值初始化为：')
            print('off_pFAK', 'on_pFAK', 'on_FAK', 'off_FAK')
            print(self.off_pFAK, ',', self.on_pFAK, ',', self.on_FAK, ',', self.off_FAK)
            # pl.plot(t, off_pFAK_list, label='off_pFAK', linewidth=3, c='#1f78b4')
            # pl.plot(t, on_pFAK_list, label='on_pFAK', linewidth=3, c='#33a02c')
            # pl.plot(t, on_FAK_list, label='on_FAK', linewidth=3, c='#e31a1c')
            # pl.plot(t, off_FAK_list, label='off_FAK', linewidth=3, c='#6a3d9a')
            # pl.legend(loc='best')
            # pl.show()

    def update(self, change=False):
        if change:
            if self.statu == 'dark':
                self.statu = 'light'
                # self.target_k2 *= self.sets_of_k[9]
                self.target_k2 = self.sets_of_k[2] * self.sets_of_k[9]
                self.change_pFAK(-self.p)
            elif self.statu == 'light':
                self.statu = 'dark'
                # self.target_k2 /= self.sets_of_k[9]
                self.target_k2 = self.sets_of_k[2] / self.sets_of_k[9]
                self.change_pFAK(self.p / (1 - self.p))

        self.all_pFAK = self.on_pFAK + self.off_pFAK

        # Way 1 off_pFAK指向on_pFAK
        way_1 = self.sets_of_k[0] * self.off_pFAK
        # Way 2 on_pFAK指向off_pFAK
        way_2 = self.sets_of_k[1] * self.on_pFAK
        # Way 3 off_FAK指向on_FAK
        way_3 = (self.sets_of_k[2] + self.sets_of_k[3] * self.all_pFAK) * self.off_FAK
        # Way 4 on_FAK指向off_FAK
        way_4 = self.sets_of_k[4] * self.on_FAK
        # Way 5 off_pFAK指向off_FAK
        way_5 = self.sets_of_k[5] * self.off_pFAK
        # Way 6 on_FAK指向on_pFAK
        way_6 = (self.sets_of_k[6] + self.sets_of_k[7] * self.all_pFAK) * self.on_FAK
        # 改变k2
        delta_k2 = self.target_k2 - self.sets_of_k[2]
        change_of_k2 = delta_k2 * self.sets_of_k[10]

        self.off_pFAK += (way_2-way_1-way_5) * self.delta_t
        self.on_pFAK += (way_1+way_6-way_2) * self.delta_t
        self.on_FAK += (way_3-way_4-way_6) * self.delta_t
        self.off_FAK += (way_4+way_5-way_3) * self.delta_t
        self.sets_of_k[2] += change_of_k2 * self.delta_t

    def get_statu(self):
        return self.off_pFAK, self.on_pFAK, self.on_FAK, self.off_FAK, self.statu

    def get_k(self):
        return self.sets_of_k

    def change_pFAK(self, n):
        self.all_pFAK = self.on_pFAK + self.off_pFAK

        delta_on_pFAK = self.on_pFAK * n
        if n <= 0:
            self.on_pFAK += delta_on_pFAK
            self.off_pFAK -= delta_on_pFAK
        else:
            part_of_off_pFAK = self.sets_of_k[0] * self.off_pFAK
            # part_of_off_pFAK = self.sets_of_k[0]
            part_of_on_FAK = (self.sets_of_k[6] + self.sets_of_k[7] * self.all_pFAK) * self.on_FAK
            # part_of_on_FAK = (self.sets_of_k[6] + self.sets_of_k[7] * self.all_pFAK)
            ratio_of_off_pFAK = part_of_off_pFAK / (part_of_off_pFAK + part_of_on_FAK)

            delta_off_pFAK = delta_on_pFAK * ratio_of_off_pFAK
            delta_on_FAK = delta_on_pFAK - delta_off_pFAK

            self.on_pFAK += delta_on_pFAK
            self.off_pFAK -= delta_off_pFAK
            self.on_FAK -= delta_on_FAK

class model_unchange:
    def __init__(self, delta_t, sets_of_k, pre_set, show=True):
        '''
        :param delta_t: 时间间隔
        :param sets_of_k:k值，单位 s-1？
        k0 off_pFAK的kin
        k1 on_pFAK的kout
        k2 off_FAK的kin，基础值
        k3 AP影响off_FAK的kin
        k4 on_FAK的kout
        k5 off_pFAK的去磷酸化
        k6 on_FAK的磷酸化，基础值
        k7 AP影响on_FAK的磷酸化
        k8 崩解率，0到0.2（0.25？）
        k9 k2的变化幅度
        k10 k2的变化速度
        :param show:展示细节
        :param preset:预设FAK量，减少计算量
        '''
        self.delta_t = delta_t
        self.sets_of_k = sets_of_k
        self.pre_set = pre_set
        self.show = show
        self.statu = 'dark'
        self.p = self.sets_of_k[8]
        self.all_pFAK = self.pre_set[0] + self.pre_set[1]
        if show:
            print('Model initializing...')
        self.get_stable()

        if show:
            print('Model Ready')
    def get_stable(self):
        if self.show:
            print('Loading...')
        self.off_pFAK = self.pre_set[0]
        self.on_pFAK = self.pre_set[1]
        self.on_FAK = self.pre_set[2]
        self.off_FAK = self.pre_set[3]
        time = 0

        t = []
        off_pFAK_list = []
        on_pFAK_list = []
        on_FAK_list = []
        off_FAK_list = []

        while True:
            if self.show:
                print('\r', time, end='', flush=True)
            off_pFAK_list.append(self.off_pFAK)
            on_pFAK_list.append(self.on_pFAK)
            on_FAK_list.append(self.on_FAK)
            off_FAK_list.append(self.off_FAK)
            t.append(time)
            self.update()
            time += self.delta_t
            if stable(off_pFAK_list, self.delta_t) and stable(on_pFAK_list, self.delta_t) \
                    and stable(on_FAK_list,self.delta_t) and stable(off_FAK_list, self.delta_t):
                break
            if min(self.off_pFAK, self.on_pFAK, self.on_FAK, self.off_FAK) < 0:
                if self.show:
                    print()
                    print('Too Low!!!')
                    print(self.off_pFAK, self.on_pFAK, self.on_FAK, self.off_FAK)
                    pl.plot(t, off_pFAK_list, label='off_pFAK', linewidth=3, c='#1f78b4')
                    pl.plot(t, on_pFAK_list, label='on_pFAK', linewidth=3, c='#33a02c')
                    pl.plot(t, on_FAK_list, label='on_FAK', linewidth=3, c='#e31a1c')
                    pl.plot(t, off_FAK_list, label='off_FAK', linewidth=3, c='#6a3d9a')
                    pl.legend(loc='best')
                    pl.show()

        if self.show:
            print()
            print('平衡态FAK值初始化为：')
            print('off_pFAK', 'on_pFAK', 'on_FAK', 'off_FAK')
            print(self.off_pFAK, ',', self.on_pFAK, ',', self.on_FAK, ',', self.off_FAK)
            # pl.plot(t, off_pFAK_list, label='off_pFAK', linewidth=3, c='#1f78b4')
            # pl.plot(t, on_pFAK_list, label='on_pFAK', linewidth=3, c='#33a02c')
            # pl.plot(t, on_FAK_list, label='on_FAK', linewidth=3, c='#e31a1c')
            # pl.plot(t, off_FAK_list, label='off_FAK', linewidth=3, c='#6a3d9a')
            # pl.legend(loc='best')
            # pl.show()
    def update(self, change=False):
        if change:
            if self.statu == 'dark':
                self.statu = 'light'
                self.change_pFAK(-self.p)
            elif self.statu == 'light':
                self.statu = 'dark'
                self.change_pFAK(self.p / (1 - self.p))

        # Way 1 off_pFAK指向on_pFAK
        way_1 = self.sets_of_k[0] * self.off_pFAK
        # Way 2 on_pFAK指向off_pFAK
        way_2 = self.sets_of_k[1] * self.on_pFAK
        # Way 3 off_FAK指向on_FAK
        way_3 = (self.sets_of_k[2] + self.sets_of_k[3] * self.all_pFAK) * self.off_FAK
        # Way 4 on_FAK指向off_FAK
        way_4 = self.sets_of_k[4] * self.on_FAK
        # Way 5 off_pFAK指向off_FAK
        way_5 = self.sets_of_k[5] * self.off_pFAK
        # Way 6 on_FAK指向on_pFAK
        way_6 = (self.sets_of_k[6] + self.sets_of_k[7] * self.all_pFAK) * self.on_FAK

        self.off_pFAK += (way_2-way_1-way_5) * self.delta_t
        self.on_pFAK += (way_1+way_6-way_2) * self.delta_t
        self.on_FAK += (way_3-way_4-way_6) * self.delta_t
        self.off_FAK += (way_4+way_5-way_3) * self.delta_t
    def get_statu(self):
        return self.off_pFAK, self.on_pFAK, self.on_FAK, self.off_FAK, self.statu

    def get_k(self):
        return self.sets_of_k

    def change_pFAK(self, n):
        delta_on_pFAK = self.on_pFAK * n
        if n <= 0:
            self.on_pFAK += delta_on_pFAK
            self.off_pFAK -= delta_on_pFAK
        else:
            part_of_off_pFAK = self.sets_of_k[0] * self.off_pFAK
            # part_of_off_pFAK = self.sets_of_k[0]
            part_of_on_FAK = (self.sets_of_k[6] + self.sets_of_k[7] * self.all_pFAK) * self.on_FAK
            # part_of_on_FAK = (self.sets_of_k[6] + self.sets_of_k[7] * self.all_pFAK)
            ratio_of_off_pFAK = part_of_off_pFAK / (part_of_off_pFAK + part_of_on_FAK)

            delta_off_pFAK = delta_on_pFAK * ratio_of_off_pFAK
            delta_on_FAK = delta_on_pFAK - delta_off_pFAK

            self.on_pFAK += delta_on_pFAK
            self.off_pFAK -= delta_off_pFAK
            self.on_FAK -= delta_on_FAK




def change_once(delta_t, sets_of_k, pre_set=0, show=False, time_limit=-1):
    A = model(delta_t, sets_of_k, show, pre_set)

    TT = []
    off_pFAK = []
    on_pFAK = []
    on_FAK = []
    off_FAK = []
    statu = []

    changed_flag = False

    t = 0

    while True:
        TT.append(t)
        new_off_pFAK, new_on_pFAK, new_on_FAK, new_off_FAK, new_statu = A.get_statu()

        off_pFAK.append(new_off_pFAK)
        on_pFAK.append(new_on_pFAK)
        on_FAK.append(new_on_FAK)
        off_FAK.append(new_off_FAK)
        statu.append(new_statu)



        if t >= delta_t*300 and changed_flag == False:
            A.update(change=True)
            changed_flag = True
        else:
            A.update(change=False)

        if time_limit > 0 and t-delta_t*300 >= time_limit:
            if show:
                print()
                print('Out of time')
            break

        if min(new_off_pFAK, new_on_pFAK, new_on_FAK, new_off_FAK) < 0:
            if show:
                print()
                print('too low')
            break
        if t>delta_t*600 and stable(off_pFAK, delta_t) and stable(on_pFAK, delta_t) and stable(on_FAK, delta_t) and stable(off_FAK, delta_t):
            if show:
                print()
                print('stable')
            break
        if show:
            print('\r', new_statu, t, end='', flush=True)
        t += delta_t
    return (TT, off_pFAK, on_pFAK, on_FAK, off_FAK, statu)

def change_cycle(delta_t, sets_of_k, cycle_time, time_limit, pre_set=0, show=False):
    A = model(delta_t, sets_of_k, show, pre_set)
    # A = model_unchange(delta_t, sets_of_k, pre_set, show)
    TT = []
    off_pFAK = []
    on_pFAK = []
    on_FAK = []
    off_FAK = []
    statu = []
    t = 0

    if time_limit > 0:  # 正数：改变几次光照
        change_time_point = np.arange(cycle_time, cycle_time * (time_limit + 0.5), cycle_time)
    else:  # 负数：总持续时间
        change_time_point = np.arange(cycle_time, (-time_limit) + 0.5 * cycle_time, cycle_time)

    while True:
        if t >= change_time_point[0]:  # 恰好达到过时间节点
            change_time_point = np.delete(change_time_point, 0, axis=0)
            if len(change_time_point) == 0:  # 结束循环？
                break
            else:  # 继续循环
                A.update(change=True)
        else:
            A.update(change=False)
        TT.append(t)
        t += delta_t

        new_off_pFAK, new_on_pFAK, new_on_FAK, new_off_FAK, new_statu = A.get_statu()
        off_pFAK.append(new_off_pFAK)
        on_pFAK.append(new_on_pFAK)
        on_FAK.append(new_on_FAK)
        off_FAK.append(new_off_FAK)
        statu.append(new_statu)
        if show:
            print('\r' + new_statu + '\t', t, end='', flush=True)
    if show:
        print()
    return TT, off_pFAK, on_pFAK, on_FAK, off_FAK, statu

def change_cycle_k(delta_t, sets_of_k, cycle_time, time_limit, pre_set=0, show=False):
    A = model(delta_t, sets_of_k, show, pre_set)

    TT = []
    off_pFAK = []
    on_pFAK = []
    on_FAK = []
    off_FAK = []
    statu = []
    k_in = []
    k_out = []
    t = 0

    if time_limit > 0:  # 正数：改变几次光照
        change_time_point = np.arange(cycle_time, cycle_time * (time_limit + 0.5), cycle_time)
    else:  # 负数：总持续时间
        change_time_point = np.arange(cycle_time, (-time_limit) + 0.5 * cycle_time, cycle_time)

    while True:
        if t >= change_time_point[0]:  # 恰好达到过时间节点
            change_time_point = np.delete(change_time_point, 0, axis=0)
            if len(change_time_point) == 0:  # 结束循环？
                break
            else:  # 继续循环
                A.update(change=True)
        else:
            A.update(change=False)
        TT.append(t)
        t += delta_t

        new_off_pFAK, new_on_pFAK, new_on_FAK, new_off_FAK, new_statu = A.get_statu()
        new_sets_of_k = A.get_k()
        if new_statu == 'dark':
            new_stiffness = 1
        else:
            new_stiffness = 0.8

        new_all_pFAK = (new_on_pFAK + new_off_pFAK)/(new_on_pFAK + new_off_pFAK + new_sets_of_k[12])
        new_k_in = new_sets_of_k[0] + new_sets_of_k[1] * new_all_pFAK * new_stiffness
        new_k_out = new_sets_of_k[2] + new_sets_of_k[3] * np.exp((new_sets_of_k[4] * new_all_pFAK* new_stiffness) / new_on_pFAK)



        off_pFAK.append(new_off_pFAK)
        on_pFAK.append(new_on_pFAK)
        on_FAK.append(new_on_FAK)
        off_FAK.append(new_off_FAK)
        statu.append(new_statu)
        k_in.append(new_k_in)
        k_out.append(new_k_out)
        if show:
            print('\r' + new_statu + '\t', t, end='', flush=True)
    if show:
        print()
    return TT, np.array(off_pFAK) + np.array(on_pFAK), np.array(k_in), np.array(k_out), statu
def stable_light_cycle(delta_t, sets_of_k, time_for_light, pre_set=0, show=False):
    # 确立模型
    A = model(delta_t, sets_of_k, show, pre_set)

    TT = []
    off_pFAK = []
    on_pFAK = []
    on_FAK = []
    off_FAK = []
    statu = []
    t = 0

    # 持续光照，时间为time_for_light
    change_flag = False
    time_of_strat = -1
    while True:
        t += delta_t
        if show:
            print('\r', t, end='', flush=True)

        TT.append(t)
        new_off_pFAK, new_on_pFAK, new_on_FAK, new_off_FAK, new_statu = A.get_statu()

        off_pFAK.append(new_off_pFAK)
        on_pFAK.append(new_on_pFAK)
        on_FAK.append(new_on_FAK)
        off_FAK.append(new_off_FAK)
        statu.append(new_statu)

        if t >= 5 and change_flag == False:
            change_flag = True
            A.update(change=True)
            time_of_strat = t
        else:
            A.update(change=False)
        if t - time_of_strat >= time_for_light and change_flag:
            if show:
                print()
                print('light over')
            break
    # 循环光照1小时
    change_time_point = np.arange(0.00, 0.6 + 0.015, 0.01)
    change_time_point += TT[-1]

    while True:
        if t >= change_time_point[0]:  # 恰好达到过时间节点
            change_time_point = np.delete(change_time_point, 0, axis=0)
            if len(change_time_point) == 0:  # 结束循环？
                break
            else:  # 继续循环
                A.update(change=True)
        else:
            A.update(change=False)
        TT.append(t)
        t += delta_t

        new_off_pFAK, new_on_pFAK, new_on_FAK, new_off_FAK, new_statu = A.get_statu()
        off_pFAK.append(new_off_pFAK)
        on_pFAK.append(new_on_pFAK)
        on_FAK.append(new_on_FAK)
        off_FAK.append(new_off_FAK)
        statu.append(new_statu)
        if show:
            print('\r' + new_statu + '\t', t, end='', flush=True)
    return TT, off_pFAK, on_pFAK, on_FAK, off_FAK, statu
def check(delta_t, sets_of_k):
    # 0 寻找平衡态
    A = model(delta_t, sets_of_k, show=False, pre_set=0)
    off_pFAK, on_pFAK, on_FAK, off_FAK, _ = A.get_statu()
    pre_set = [off_pFAK, on_pFAK, on_FAK, off_FAK]

    # 检验1 黑暗培养平衡后循环光照固定时间
    result_1 = []
    # 1）持续黑暗（跳过，作为基准当然是1）       1
    pass
    # 2）持续光照至平滑          0.825977718976908
    TT, off_pFAK, on_pFAK, on_FAK, off_FAK = change_once(delta_t, sets_of_k, pre_set=0, show=False, time_limit=-1)
    result_1.append(((off_pFAK[-1] + on_pFAK[-1]) / (off_pFAK[0] + on_pFAK[0])))
    # 3）1min间隔循环持续12h（720min）    3.08483936663155
    TT, off_pFAK, on_pFAK, on_FAK, off_FAK = change_cycle(delta_t, sets_of_k, cycle_time=60, time_limit=-720 * 60,
                                                          pre_set=0, show=False)
    result_1.append(((off_pFAK[-1] + on_pFAK[-1]) / (off_pFAK[0] + on_pFAK[0])))
    # 4）5min 同上                  1.65986761872456
    TT, off_pFAK, on_pFAK, on_FAK, off_FAK = change_cycle(delta_t, sets_of_k, cycle_time=300, time_limit=-720 * 60,
                                                          pre_set=0, show=False)
    result_1.append(((off_pFAK[-1] + on_pFAK[-1]) / (off_pFAK[0] + on_pFAK[0])))
    # 5）30min 同上                 0.96247625260154
    TT, off_pFAK, on_pFAK, on_FAK, off_FAK = change_cycle(delta_t, sets_of_k, cycle_time=1800, time_limit=-720 * 60,
                                                          pre_set=0, show=False)
    result_1.append(((off_pFAK[-1] + on_pFAK[-1]) / (off_pFAK[0] + on_pFAK[0])))

    # 检验2 黑暗培养平衡后循环光照固定次数(96次）
    result_2 = []
    # 1）1min 间隔      1.79340534666633
    TT, off_pFAK, on_pFAK, on_FAK, off_FAK = change_cycle(delta_t, sets_of_k, cycle_time=60, time_limit=96, pre_set=0,
                                                          show=False)
    result_2.append(((off_pFAK[-1] + on_pFAK[-1]) / (off_pFAK[0] + on_pFAK[0])))
    # 2）5min 间隔      2.17397684802067
    TT, off_pFAK, on_pFAK, on_FAK, off_FAK = change_cycle(delta_t, sets_of_k, cycle_time=300, time_limit=96, pre_set=0,
                                                          show=False)
    result_2.append(((off_pFAK[-1] + on_pFAK[-1]) / (off_pFAK[0] + on_pFAK[0])))
    # 1）30min 间隔     1
    TT, off_pFAK, on_pFAK, on_FAK, off_FAK = change_cycle(delta_t, sets_of_k, cycle_time=1800, time_limit=96, pre_set=0,
                                                          show=False)
    result_2.append(((off_pFAK[-1] + on_pFAK[-1]) / (off_pFAK[0] + on_pFAK[0])))

    # 检验3 黑暗培养平衡后1min循环光照
    result_3 = []
    TT, off_pFAK, on_pFAK, on_FAK, off_FAK = change_cycle(delta_t, sets_of_k, cycle_time=60, time_limit=-720 * 60,
                                                          pre_set=0,
                                                          show=False)
    # 1) 1h      1.41894159015926
    result_3.append(
        ((off_pFAK[-round(660 * 60 / delta_t)] + on_pFAK[-round(660 / delta_t)]) / (off_pFAK[0] + on_pFAK[0])))
    # 2) 3h      2.17691536267287
    result_3.append(
        ((off_pFAK[-round(540 * 60 / delta_t)] + on_pFAK[-round(540 / delta_t)]) / (off_pFAK[0] + on_pFAK[0])))
    # 3) 6h      3.39073421878973
    result_3.append(
        ((off_pFAK[-round(360 * 60 / delta_t)] + on_pFAK[-round(360 / delta_t)]) / (off_pFAK[0] + on_pFAK[0])))
    # 4) 12h     4.07846763136423
    result_3.append(((off_pFAK[-1] + on_pFAK[-1]) / (off_pFAK[0] + on_pFAK[0])))

    # 检验4 黑暗培养后持续光照
    result_4 = []
    TT, off_pFAK, on_pFAK, on_FAK, off_FAK = change_once(delta_t, sets_of_k, pre_set=0, show=False, time_limit=-1)
    # 1) 持续光照1min    0.965441557810043
    result_4.append(((off_pFAK[round(60 + 5 / delta_t)] + on_pFAK[round(60 / delta_t)]) / (off_pFAK[0] + on_pFAK[0])))
    # 2) 持续光照5min    0.876029243828745
    result_4.append(
        ((off_pFAK[round(300 + 5 / delta_t)] + on_pFAK[round(3000 / delta_t)]) / (off_pFAK[0] + on_pFAK[0])))
    # 1) 持续光照30min   0.630281601532441
    result_4.append(
        ((off_pFAK[round(1800 + 5 / delta_t)] + on_pFAK[round(1800 + 5 / delta_t)]) / (off_pFAK[0] + on_pFAK[0])))

    # 检验5 先黑岩平衡，再持续光照不同时间，最后1min频率循环1h
    result_5 = []
    # 1) 光照时间1min    2.45860150474854
    TT, off_pFAK, on_pFAK, on_FAK, off_FAK = stable_light_cycle(delta_t, sets_of_k, 60, pre_set=0, show=False)
    result_5.append(((off_pFAK[-1] + on_pFAK[-1]) / (off_pFAK[0] + on_pFAK[0])))

    # 2) 光照时间5min    2.35842271464659
    TT, off_pFAK, on_pFAK, on_FAK, off_FAK = stable_light_cycle(delta_t, sets_of_k, 300, pre_set=0, show=False)
    result_5.append(((off_pFAK[-1] + on_pFAK[-1]) / (off_pFAK[0] + on_pFAK[0])))

    # 3) 光照时间30min   1.77304244832633
    TT, off_pFAK, on_pFAK, on_FAK, off_FAK = stable_light_cycle(delta_t, sets_of_k, 1800, pre_set=0, show=False)
    result_5.append(((off_pFAK[-1] + on_pFAK[-1]) / (off_pFAK[0] + on_pFAK[0])))

    return result_1, result_2, result_3, result_4, result_5
def estimate_k(sets_of_k, show=False):
    A = model(0.1, np.copy(sets_of_k), show=False, pre_set=0)
    new_off_pFAK, new_on_pFAK, new_on_FAK, new_off_FAK, _ = A.get_statu()
    pre_set = [new_off_pFAK, new_on_pFAK, new_on_FAK, new_off_FAK]

    result_1 = []
    TT, off_pFAK, on_pFAK, on_FAK, off_FAK, statu = change_cycle(0.1, np.copy(sets_of_k), cycle_time=1,
                                                                 time_limit=720, pre_set=pre_set, show=False)
    all_pFAK = np.array(off_pFAK) + np.array(on_pFAK)
    all_pFAK = all_pFAK / all_pFAK[0]
    result_1.append(all_pFAK[600-1])
    result_1.append(all_pFAK[1800-1])
    result_1.append(all_pFAK[3600-1])
    result_1.append(all_pFAK[7200-1])
    result_2 = [result_1[-1]]
    for cycle_time in [5,30]:
        TT, off_pFAK, on_pFAK, on_FAK, off_FAK,statu = change_cycle(0.1, np.copy(sets_of_k), cycle_time=cycle_time, time_limit=-720, pre_set=pre_set, show=False)
        result_2.append((off_pFAK[-1]+on_pFAK[-1]) / (off_pFAK[0]+on_pFAK[0]))
    if show:
        print(result_1)
        print(result_2)
    return np.array(result_1),np.array(result_2)
def mark_K(sets_of_k, show=False):
    result_1, result_2 = estimate_k(sets_of_k,show)
    gt_1 = np.array([1.41894159, 2.176915363, 3.390734219, 4.078467631])
    gt_2 = np.array([3.084839367, 1.659867619, 0.962476253])/0.962476253
    error_1 = np.power(result_1-gt_1, 2)
    error_2 = np.power(result_2-gt_2, 2)
    # error[1] *= 3
    error_r = np.sum(error_1) + np.sum(error_2)
    # error_r += ((result[0]-result[1])*3)
    return error_r
if __name__ == '__main__':
    # a = open('stiffness_pFAK.txt', 'w', encoding='utf8')
    # factor_range = np.arange(0.5,2.01,0.1)
    # yy = []
    # for factor in factor_range:
    #     print(factor)
    #     a.write(str(factor)+'\t')
    #     k0 = 0.163 * 0.87  # off_pFAK的kin 0.14181
    #     k1 = 0.0834 * 1.13  # on_pFAK的kout    0.094242
    #     k2_b = 0.030858525
    #     k2_s = 0.038486475*factor
    #     k2 = k2_b + k2_s
    #     # k2 = 0.069345  # off_FAK的kin，基础值
    #     k3 = 0.000012495  # AP影响off_FAK的kin
    #     k4 = 0.0958  # on_FAK的kout
    #     k5 = 0.07917648  # off_pFAK的去磷酸化
    #     k6 = 0.0000649  # on_FAK的磷酸化，基础值        new k0
    #     k7 = 0.0000570  # AP影响on_FAK的磷酸化
    #     k8 = 0.180  # 崩解率，0到0.2（0.25？）
    #     k9 = (k2_b + 0.8 * k2_s) / k2
    #     # k9 = 0.889  # k2的变化幅度
    #     k10 = 0.138  # k2的变化速度
    #
    #     sets_of_k = np.array([k0, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10])
    #     A = model(0.1, np.copy(sets_of_k), show=False, pre_set=0)
    #     new_off_pFAK, new_on_pFAK, new_on_FAK, new_off_FAK, _ = A.get_statu()
    #     pre_set = [new_off_pFAK, new_on_pFAK, new_on_FAK, new_off_FAK]
    #     TT, off_pFAK, on_pFAK, on_FAK, off_FAK, statu = change_cycle(0.01, np.copy(sets_of_k), cycle_time=1,
    #                                                                  time_limit=720, pre_set=pre_set, show=False)
    #     all_pFAK = np.array(off_pFAK) + np.array(on_pFAK)
    #     y=all_pFAK[-1]
    #     yy.append(y)
    #     a.write(str(y)+'\n')
    # a.close()
    # pl.plot(factor_range,yy)
    # pl.show()





    # delta_t = 0.0001
    #
    k0 = 0.163*0.87  # off_pFAK的kin 0.14181
    # k0 *= 0.5
    # k0 = 0.163*0.77  # off_pFAK的kin 0.14181
    k1 = 0.0834*1.13  # on_pFAK的kout    0.094242
    # k1 *= 0.5
    # k1 = 0.0834*1.03  # on_pFAK的kout    0.094242

    k2_b = 0.030858525
    k2_s = 0.038486475
    k2 = k2_b + k2_s
    # k2 *= 0.5
    # k2 = 0.069345  # off_FAK的kin，基础值

    k3 = 0.000012495 # AP影响off_FAK的kin
    # k3 *= 0.5
    k4 = 0.0958 # on_FAK的kout
    k5 = 0.07917648  # off_pFAK的去磷酸化
    k6 = 0.0000649 # on_FAK的磷酸化，基础值        new k0
    # k6 *= 0.5
    k7 = 0.0000570  # AP影响on_FAK的磷酸化
    # k7 *= 0.5
    k8 = 0.180  # 崩解率，0到0.2（0.25？）

    k9 = (k2_b + 0.8*k2_s)/k2
    # k9 = 0.889  # k2的变化幅度

    k10 = 0.138  # k2的变化速度

    sets_of_k = np.array([k0, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10])
    # sets_of_k = np.array([0.16268250974369255, 0.08551437822805749, 0.04972087284801394, 7.09561432848036e-06, 0.09874930968536277, 0.08189682471261864, 6.441979967161009e-05, 5.80448949616922e-05, 0.2, 1.0, 0.0])
    # f = open('Dephosphorylation.txt', 'w', encoding='utf8')
    # for factor in np.arange(0.1,10.01,0.1):
    #     print(factor)
    #     test_k = np.copy(sets_of_k)
    #     test_k[5] *= factor
    #     A = model(0.1, np.copy(sets_of_k), show=False, pre_set=0)
    #     new_off_pFAK, new_on_pFAK, new_on_FAK, new_off_FAK, _ = A.get_statu()
    #     pre_set = [new_off_pFAK, new_on_pFAK, new_on_FAK, new_off_FAK]
    #     TT, off_pFAK, on_pFAK, on_FAK, off_FAK, statu = change_cycle(0.1, np.copy(test_k), cycle_time=1,
    #                                                                  time_limit=720, pre_set=pre_set, show=False)
    #     all_pFAK = np.array(off_pFAK) + np.array(on_pFAK)
    #     f.write(str(factor) + '\t' + str(all_pFAK[-1]/all_pFAK[0]) + '\n')
    # f.close()



    # pre_set = [75.1671813385021, 104.31826641130729, 89.51215557027695, 731.0023966799039]
    #
    A = model(0.1, np.copy(sets_of_k), show=False, pre_set=0)
    new_off_pFAK, new_on_pFAK, new_on_FAK, new_off_FAK, _ = A.get_statu()
    pre_set = [new_off_pFAK, new_on_pFAK, new_on_FAK, new_off_FAK]
    # print (pre_set)

    #1min 720cycles

    TT, off_pFAK, on_pFAK, on_FAK, off_FAK, statu = change_cycle(0.01, np.copy(sets_of_k), cycle_time=1,
                                                                     time_limit=10000, pre_set=pre_set, show=False)
    all_pFAK = np.array(off_pFAK) + np.array(on_pFAK)
    all_FAK = np.array(on_FAK) + np.array(off_FAK)
    # all_pFAK = all_pFAK/all_pFAK[0]
    # x = len(all_pFAK)
    # f = open('1min_for_20m_nochanged.txt', 'w', encoding='utf8')
    # for i in range(x):
    #     f.write(str(TT[i]) + '\t' + str(all_pFAK[i]) + '\n')
    # for i in np.arange(0,721,2):
    #     print(i)
    #     y = all_pFAK[max(i*100-1, 0)]
    #     f.write(str(i)+'\t'+str(y)+'\n')
    # f.close()
    pl.plot(TT, all_pFAK, label='pFAK')
    pl.plot(TT, all_FAK, label='FAK')
    pl.legend(loc='best')
    pl.show()

    # print (len(all_pFAK))
    # print("on_pFAK:",str(np.array(on_pFAK)[-1]/np.array(on_pFAK)[0]))
    # print (all_pFAK[6000-1])
    # print (all_pFAK[18000-1])
    # print (all_pFAK[36000-1])
    # print (all_pFAK[72000-1])
    # print('---------------------')
    # print(all_pFAK[-1] / all_pFAK[0])
    # TT, off_pFAK, on_pFAK, _, _, _ = change_cycle(0.01, np.copy(sets_of_k), cycle_time=1,
    #                                                                  time_limit=96, pre_set=pre_set, show=False)
    # all_pFAK = np.array(off_pFAK) + np.array(on_pFAK)
    # print(all_pFAK[-1] / all_pFAK[0])
    # TT, off_pFAK, on_pFAK, _, _, _ = change_cycle(0.1, np.copy(sets_of_k), cycle_time=5,
    #                                                                  time_limit=96, pre_set=pre_set, show=False)
    # all_pFAK = np.array(off_pFAK) + np.array(on_pFAK)
    # print(all_pFAK[-1] / all_pFAK[0])
    # TT, off_pFAK, on_pFAK, _, _, _ = change_cycle(0.1, np.copy(sets_of_k), cycle_time=30,
    #                                                                  time_limit=96, pre_set=pre_set, show=False)
    # all_pFAK = np.array(off_pFAK) + np.array(on_pFAK)
    # print(all_pFAK[-1] / all_pFAK[0])
    # print('---------------------')
    '''
    0.1到30分钟，96次循环
    
    t_range = np.arange(0.1,30.1,0.1)
    f = open('96_cycles.txt', 'w', encoding='utf8')
    for t in t_range:
        print(t)
        if t <= 1:
            cycle_time = 0.001
        elif t <= 10:
            cycle_time = 0.01
        else:
            cycle_time = 0.1
        TT, off_pFAK, on_pFAK, on_FAK, off_FAK, statu = change_cycle(cycle_time, np.copy(sets_of_k), cycle_time=t,
                                                                     time_limit=96, pre_set=pre_set, show=False)
        all_pFAK = np.array(off_pFAK) + np.array(on_pFAK)
        f.write(str(t) + '\t' + str(all_pFAK[-1]/all_pFAK[0]) + '\n')
    f.close()
    '''
    # for t in t_range_2:
    #     # print(t)
    #     TT, off_pFAK, on_pFAK, on_FAK, off_FAK, statu = change_cycle(0.1, np.copy(sets_of_k), cycle_time=t,
    #                                                                  time_limit=96, pre_set=pre_set, show=False)
    #     all_pFAK = np.array(off_pFAK) + np.array(on_pFAK)
    #     result.append(all_pFAK[-1]/all_pFAK[0])
    #     tt.append(t)
    # pl.plot(tt, result)
    # pl.show()
    #
    # # result = []
    # # t_range = np.arange(0.1, 1, 0.1)
    # # for t in t_range:
    # #     print(t)
    # #     TT, off_pFAK, on_pFAK, on_FAK, off_FAK,statu = change_cycle(0.001, np.copy(sets_of_k), cycle_time=t, time_limit=96, pre_set=pre_set, show=False)
    # #     result.append ((off_pFAK[-1]+on_pFAK[-1])/(off_pFAK[0]+on_pFAK[0]))
    # # pl.plot(t_range, result)
    # # pl.show()
    #
    # # TT, all_pFAK, k_in, k_out, statu = change_cycle_k(0.1, np.copy(sets_of_k), 1, 720, pre_set=pre_set, show=False)
    # # pl.plot(TT,all_pFAK,label='all_pFAK')
    # # pl.show()
    # # pl.plot(TT,k_in,label='in')
    # # pl.show()
    # # pl.plot(TT,k_out,label='out')
    # # pl.show()
    # # #
    # # pl.plot(TT,k_in/k_in[0],label='in')
    # # pl.plot(TT, all_pFAK/all_pFAK[0], label='allpFAK')
    # # pl.plot(TT,k_out/k_out[0],label='out')
    # # pl.legend(loc='best')
    # # pl.show()
    #
    # # TT, off_pFAK, on_pFAK, on_FAK, off_FAK, statu = change_once(0.1, sets_of_k, pre_set=pre_set, show=True, time_limit=-1)
    # # TT, off_pFAK, on_pFAK, on_FAK, off_FAK, statu = stable_light_cycle(delta_t, sets_of_k, 0.3, pre_set=pre_set,show=True)
    # # TT, off_pFAK, on_pFAK, on_FAK, off_FAK, statu = change_once(0.1, np.copy(sets_of_k), pre_set=pre_set, show=True, time_limit=-1)
    # # pl.plot(TT, (np.array(on_pFAK) + np.array(off_pFAK)), label='all_pFAK', linewidth=3, c='#b15928')
    # # pl.fill_between(TT, np.zeros_like(TT)+np.max(np.array(on_pFAK) + np.array(off_pFAK)), np.zeros_like(TT)+np.min(np.array(on_pFAK) + np.array(off_pFAK)), where=np.array(statu)=='dark',alpha=0.2)
    # # pl.legend(loc='best')
    # # all_pFAK = np.array(off_pFAK) + np.array(on_pFAK)
    # # print()
    # # print(all_pFAK[-1] / all_pFAK[0])
    # # pl.show()
    #
    '''
    1到360分钟，720分钟循环
    
    t_range = np.arange(1,361,1)
    f = open('cycle_of_720m.txt', 'w', encoding='utf8')
    for t in t_range:
        if 360%t != 0:
            continue
        print(t)
        TT, off_pFAK, on_pFAK, on_FAK, off_FAK, statu = change_cycle(0.1, np.copy(sets_of_k), cycle_time=t,
                                                                     time_limit=-720, pre_set=pre_set, show=False)
        all_pFAK = np.array(off_pFAK) + np.array(on_pFAK)
        f.write(str(t) + '\t' + str(all_pFAK[-1]/all_pFAK[0]) + '\n')
    f.close()
    '''
    '''
    改变stiffness0.75 1 1.5倍率，1min 720次循环,记录allpFAK和onpFAK
    
    k0 = 0.163*0.87  # off_pFAK的kin
    k1 = 0.0834*1.13  # on_pFAK的kout
    stiff = 0.038486475
    k_base= 0.030858525
    k2 = stiff+k_base
    k3 = 0.000012495 # AP影响off_FAK的kin
    k4 = 0.0958 # on_FAK的kout
    k5 = 0.07917648  # off_pFAK的去磷酸化
    k6 = 0.0000649 # on_FAK的磷酸化，基础值
    k7 = 0.0000570  # AP影响on_FAK的磷酸化
    k8 = 0.180  # 崩解率，0到0.2（0.25？）
    k9 = (stiff*0.8+k_base)/(stiff+k_base)
    k10 = 0.138  # k2的变化速度
    sets_of_k = np.array([k0, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10])

    A = model(0.1, np.copy(sets_of_k), show=False, pre_set=0)
    new_off_pFAK, new_on_pFAK, new_on_FAK, new_off_FAK, _ = A.get_statu()
    pre_set = [new_off_pFAK, new_on_pFAK, new_on_FAK, new_off_FAK]
    TT, off_pFAK, on_pFAK, on_FAK, off_FAK, statu = change_cycle(0.01, np.copy(sets_of_k), cycle_time=1,
                                                                     time_limit=720, pre_set=pre_set, show=False)

    all_pFAK = np.array(off_pFAK) + np.array(on_pFAK)
    all_pFAK = all_pFAK/all_pFAK[0]
    on_pFAK = np.array(on_pFAK)/on_pFAK[0]
    print(len(on_pFAK))
    f = open('1 stiffness.txt', 'w', encoding='utf8')
    f.write('Time\tall_pFAk\ton_pFAK\n')
    for i in np.arange(0,721,2):
        print(i)
        y1 = all_pFAK[max(i*100-1, 0)]
        y2 = on_pFAK[max(i*100-1, 0)]
        f.write(str(i)+'\t'+str(y1)+'\t'+str(y2)+'\n')
    f.close()

    k0 = 0.163*0.87  # off_pFAK的kin
    k1 = 0.0834*1.13  # on_pFAK的kout
    stiff = 0.038486475
    stiff *= 0.75
    k_base= 0.030858525
    k2 = stiff+k_base
    k3 = 0.000012495 # AP影响off_FAK的kin
    k4 = 0.0958 # on_FAK的kout
    k5 = 0.07917648  # off_pFAK的去磷酸化
    k6 = 0.0000649 # on_FAK的磷酸化，基础值
    k7 = 0.0000570  # AP影响on_FAK的磷酸化
    k8 = 0.180  # 崩解率，0到0.2（0.25？）
    k9 = (stiff*0.8+k_base)/(stiff+k_base)
    k10 = 0.138  # k2的变化速度
    sets_of_k = np.array([k0, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10])

    A = model(0.1, np.copy(sets_of_k), show=False, pre_set=0)
    new_off_pFAK, new_on_pFAK, new_on_FAK, new_off_FAK, _ = A.get_statu()
    pre_set = [new_off_pFAK, new_on_pFAK, new_on_FAK, new_off_FAK]
    TT, off_pFAK, on_pFAK, on_FAK, off_FAK, statu = change_cycle(0.01, np.copy(sets_of_k), cycle_time=1,
                                                                     time_limit=720, pre_set=pre_set, show=False)

    all_pFAK = np.array(off_pFAK) + np.array(on_pFAK)
    all_pFAK = all_pFAK/all_pFAK[0]
    on_pFAK = np.array(on_pFAK)/on_pFAK[0]
    f = open('0.75 stiffness.txt', 'w', encoding='utf8')
    f.write('Time\tall_pFAk\ton_pFAK\n')
    for i in np.arange(0,721,2):
        print(i)
        y1 = all_pFAK[max(i*100-1, 0)]
        y2 = on_pFAK[max(i*100-1, 0)]
        f.write(str(i)+'\t'+str(y1)+'\t'+str(y2)+'\n')
    f.close()

    k0 = 0.163*0.87  # off_pFAK的kin
    k1 = 0.0834*1.13  # on_pFAK的kout
    stiff = 0.038486475
    stiff *= 1.5
    k_base= 0.030858525
    k2 = stiff+k_base
    k3 = 0.000012495 # AP影响off_FAK的kin
    k4 = 0.0958 # on_FAK的kout
    k5 = 0.07917648  # off_pFAK的去磷酸化
    k6 = 0.0000649 # on_FAK的磷酸化，基础值
    k7 = 0.0000570  # AP影响on_FAK的磷酸化
    k8 = 0.180  # 崩解率，0到0.2（0.25？）
    k9 = (stiff*0.8+k_base)/(stiff+k_base)
    k10 = 0.138  # k2的变化速度
    sets_of_k = np.array([k0, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10])

    A = model(0.1, np.copy(sets_of_k), show=False, pre_set=0)
    new_off_pFAK, new_on_pFAK, new_on_FAK, new_off_FAK, _ = A.get_statu()
    pre_set = [new_off_pFAK, new_on_pFAK, new_on_FAK, new_off_FAK]
    TT, off_pFAK, on_pFAK, on_FAK, off_FAK, statu = change_cycle(0.01, np.copy(sets_of_k), cycle_time=1,
                                                                     time_limit=720, pre_set=pre_set, show=False)
    all_pFAK = np.array(off_pFAK) + np.array(on_pFAK)
    all_pFAK = all_pFAK/all_pFAK[0]
    on_pFAK = np.array(on_pFAK)/on_pFAK[0]
    f = open('1.5 stiffness.txt', 'w', encoding='utf8')
    f.write('Time\tall_pFAk\ton_pFAK\n')
    for i in np.arange(0,721,2):
        print(i)
        y1 = all_pFAK[max(i*100-1, 0)]
        y2 = on_pFAK[max(i*100-1, 0)]
        f.write(str(i)+'\t'+str(y1)+'\t'+str(y2)+'\n')
    f.close()
    '''

    # TT, off_pFAK, on_pFAK, on_FAK, off_FAK,statu = change_cycle(0.1, np.copy(sets_of_k), cycle_time=1, time_limit=-720, pre_set=pre_set, show=False)
    # TT = np.array(TT)/60
    # all_pFAK = np.array(off_pFAK) + np.array(on_pFAK)
    # pl.plot(TT, all_pFAK/all_pFAK[0])
    # pl.fill_between(TT, np.zeros_like(TT) + np.max(all_pFAK/all_pFAK[0]),
    #                 np.zeros_like(TT) + np.min(all_pFAK/all_pFAK[0]), where=np.array(statu) == 'dark',
    #                 alpha=0.2)
    # # print()
    # all_pFAK = all_pFAK/all_pFAK[0]
    # print(all_pFAK[-1])
    # pl.show()
    # TT, off_pFAK, on_pFAK, on_FAK, off_FAK,statu = change_cycle(0.1, np.copy(sets_of_k), cycle_time=5, time_limit=-720, pre_set=pre_set, show=False)
    # all_pFAK = np.array(off_pFAK) + np.array(on_pFAK)
    # pl.plot(TT, all_pFAK)
    # pl.fill_between(TT, np.zeros_like(TT) + np.max(all_pFAK),
    #                 np.zeros_like(TT) + np.min(all_pFAK), where=np.array(statu) == 'dark',
    #                 alpha=0.2)
    # # print()
    # print(all_pFAK[-1] / all_pFAK[0])
    # pl.show()
    # TT, off_pFAK, on_pFAK, on_FAK, off_FAK,statu = change_cycle(0.1, np.copy(sets_of_k), cycle_time=30, time_limit=-720, pre_set=pre_set, show=False)
    # all_pFAK = np.array(off_pFAK) + np.array(on_pFAK)
    # pl.plot(TT, all_pFAK)
    # pl.fill_between(TT, np.zeros_like(TT) + np.max(all_pFAK),
    #                 np.zeros_like(TT) + np.min(all_pFAK), where=np.array(statu) == 'dark',
    #                 alpha=0.2)
    # # print()
    # print(all_pFAK[-1] / all_pFAK[0])
    # pl.show()
    #
    # print (estimate_k(np.copy(sets_of_k)))


    # TT, off_pFAK, on_pFAK, on_FAK, off_FAK,statu = change_cycle(0.1, np.copy(sets_of_k), cycle_time=1, time_limit=3, pre_set=pre_set, show=True)
    # all_pFAK = np.array(off_pFAK) + np.array(on_pFAK)
    # pl.plot(TT, all_pFAK)
    # pl.fill_between(TT, np.zeros_like(TT) + np.max(all_pFAK),
    #                 np.zeros_like(TT) + np.min(all_pFAK), where=np.array(statu) == 'dark',
    #                 alpha=0.2)
    # print()
    # print(all_pFAK[-1] / all_pFAK[0])
    # pl.show()
    # TT, off_pFAK, on_pFAK, on_FAK, off_FAK,statu = change_cycle(0.01, np.copy(sets_of_k), cycle_time=0.5, time_limit=720, pre_set=pre_set, show=True)
    # all_pFAK = np.array(off_pFAK) + np.array(on_pFAK)
    # pl.plot(TT, all_pFAK)
    # pl.fill_between(TT, np.zeros_like(TT) + np.max(all_pFAK),
    #                 np.zeros_like(TT) + np.min(all_pFAK), where=np.array(statu) == 'dark',
    #                 alpha=0.2)
    # print()
    # print(all_pFAK[-1] / all_pFAK[0])
    # pl.show()
    #
    # pl.plot(TT, off_pFAK, label='off_pFAK', linewidth=3, c='#1f78b4')
    # pl.plot(TT, on_pFAK, label='on_pFAK', linewidth=3, c='#33a02c')
    # pl.plot(TT, on_FAK, label='on_FAK', linewidth=3, c='#e31a1c')
    # pl.plot(TT, off_FAK, label='off_FAK', linewidth=3, c='#6a3d9a')
    # pl.legend(loc='best')
    # top = np.max([np.max(off_pFAK),np.max(on_pFAK),np.max(on_FAK),np.max(off_FAK)])
    # bottom = np.min([np.min(off_pFAK),np.min(on_pFAK),np.min(on_FAK),np.min(off_FAK)])
    # pl.fill_between(TT, np.zeros_like(TT) + top,
    #                 np.zeros_like(TT) + bottom, where=np.array(statu) == 'dark',
    #                 alpha=0.2)
    # pl.show()
    # length_of_1m = round(0.01 /delta_t)
    # length_of_front = round(5/delta_t)
    # length_of_1h = length_of_1m * 60
    # print(length_of_1h)
    # print(len(all_pFAK))

    # print(all_pFAK[length_of_front+1*length_of_1m]/all_pFAK[0])
    # print(all_pFAK[length_of_front+5*length_of_1m]/all_pFAK[0])
    # print(all_pFAK[length_of_front+30*length_of_1m]/all_pFAK[0])


    # #
    # pl.plot(TT,k_in,label='in')
    # pl.show()
    # pl.plot(TT,k_out,label='out')
    # pl.show()
    # #
    # pl.plot(TT,k_in/k_in[0],label='in')
    # pl.plot(TT, all_pFAK/all_pFAK[0], label='allpFAK')
    # pl.plot(TT,k_out/k_out[0],label='out')
    # pl.legend(loc='best')
    # pl.show()
    #
    # length = 200
    # TT_a = np.array(TT[-72000:]).reshape(-1,length)[:,0]
    # all_pFAK_a = all_pFAK[-72000:].reshape(-1,length)[:,0]
    # k_in_a = np.array(k_in[-72000:]).reshape(-1,length)[:,0]
    # k_out_a = np.array(k_out[-72000:]).reshape(-1,length)[:,0]
    # pl.plot(TT_a,k_in_a/k_in_a[0],label='in')
    # pl.plot(TT_a, all_pFAK_a/all_pFAK_a[0], label='allpFAK')
    # pl.plot(TT_a,k_out_a/k_out_a[0],label='out')
    # pl.legend(loc='best')
    # pl.show()
    #

    # # cycle_time_range = [0.01, 0.05,0.1,0.2,0.3,0.4] + list(np.arange(0.5,50.5,0.5)) + list(np.arange(50, 500, 5))
    # pl.plot(TT, (np.array(on_pFAK) + np.array(off_pFAK)), label='all_pFAK', linewidth=3, c='#b15928')
    # pl.fill_between(TT, np.zeros_like(TT) + np.max(np.array(on_pFAK) + np.array(off_pFAK)),
    #                 np.zeros_like(TT) + np.min(np.array(on_pFAK) + np.array(off_pFAK)), where=np.array(statu) == 'dark',
    #                 alpha=0.2)
    # pl.legend(loc='best')
    # pl.show()
    #
    # A = model(delta_t*10, np.copy(sets_of_k), show=True, pre_set=0)
    # new_off_pFAK, new_on_pFAK, new_on_FAK, new_off_FAK, _ = A.get_statu()
    # pre_set = [new_off_pFAK, new_on_pFAK, new_on_FAK, new_off_FAK]
    # print (pre_set)
    #
    # TT, off_pFAK, on_pFAK, on_FAK, off_FAK, statu = change_once(delta_t*10, np.copy(sets_of_k), pre_set=pre_set, show=True, time_limit=-1)
    # pl.plot(TT, (np.array(on_pFAK) + np.array(off_pFAK)), label='all_pFAK', linewidth=3, c='#b15928')
    # pl.fill_between(TT, np.zeros_like(TT)+np.max(np.array(on_pFAK) + np.array(off_pFAK)), np.zeros_like(TT)+np.min(np.array(on_pFAK) + np.array(off_pFAK)), where=np.array(statu)=='dark',alpha=0.2)
    # pl.legend(loc='best')
    # pl.show()
    # #
    # cycle_time_range = np.arange(1,31,1)
    # result = []
    # for time in cycle_time_range:
    #     print(time)
    #     TT, off_pFAK, on_pFAK, on_FAK, off_FAK, statu = change_cycle(delta_t, np.copy(sets_of_k), cycle_time=time, time_limit=96, pre_set=pre_set, show=False)
    #     result.append ((off_pFAK[-1]+on_pFAK[-1])/(off_pFAK[0]+on_pFAK[0]))
    # pl.cla()
    # pl.plot(cycle_time_range, result)
    # pl.show()

    #
    # pl.plot(TT, off_pFAK, label='off_pFAK', linewidth=3, c='#1f78b4')
    # pl.plot(TT, on_pFAK, label='on_pFAK', linewidth=3, c='#33a02c')
    # pl.plot(TT, on_FAK, label='on_FAK', linewidth=3, c='#e31a1c')
    # pl.plot(TT, off_FAK, label='off_FAK', linewidth=3, c='#6a3d9a')
    # print ((on_pFAK[-1]+ off_pFAK[-1])/(on_pFAK[0]+ off_pFAK[0]))
    # pl.plot(TT, (np.array(on_pFAK) + np.array(off_pFAK)), label='all_pFAK', linewidth=3, c='#b15928')
    # # pl.fill_between(TT, np.zeros_like(TT) + np.max(np.array(on_pFAK) + np.array(off_pFAK)),
    # #                 np.zeros_like(TT) + np.min(np.array(on_pFAK) + np.array(off_pFAK)), where=np.array(statu) == 'dark',
    # #                 alpha=0.2)
    # pl.legend(loc='best')
    # pl.show()


    # lr = 0.00000005
    # train_time = 0
    # name = time.strftime("%Y_%m_%d_%H_%M_%S", time.localtime())
    # sets_of_k_stiff = np.copy(sets_of_k)
    # sets_of_k_stiff = np.array([0.16268250974369255, 0.08551437822805749, 0.04972087284801394, 7.09561432848036e-06, 0.09874930968536277, 0.08189682471261864, 6.441979967161009e-05, 5.80448949616922e-05, 0.2, 1.0, 0.0])
    # error_list = []
    # best_k = []
    # while True:
    #     print(train_time, sets_of_k_stiff)
    #     base_error = mark_K(sets_of_k_stiff, show=True)
    #     print(base_error, lr)
    #     if train_time == 0:
    #         mode = 'w'
    #         error_list.append(base_error)
    #         best_k = sets_of_k_stiff
    #     else:
    #         mode = 'a'
    #         if base_error < min(error_list):
    #             error_list = [base_error]
    #             best_k = sets_of_k_stiff
    #         else:
    #             error_list.append(base_error)
    #             if len(error_list) > 50:
    #                 error_list = error_list[0]
    #                 sets_of_k_stiff = best_k
    #                 base_error = mark_K(sets_of_k_stiff)
    #                 lr = lr/10
    #     f = open('result_' + name + '.txt', mode, encoding='utf8')
    #     f.write(str(train_time) + '\t' + str(lr) + '\n')
    #     f.write(str(best_k) +'\t' + str(error_list[0]) + '\n')
    #     f.write(str(base_error) + '\n')
    #     f.write(str(sets_of_k_stiff) + '\n')
    #     f.write('----------------------------\n')
    #     f.close()
    #     new_k = []
    #     for i in range(8):
    #         change = min(1e-12, lr/1000)
    #         test_k = np.copy(sets_of_k_stiff)
    #         test_k[i] = test_k[i]+change
    #         test_error = mark_K(test_k)
    #         delta_error = (test_error-base_error)/change
    #         delta_k = delta_error*lr*sets_of_k_stiff[i]
    #         # if i == 6 or i == 3 or i == 7:
    #         #     delta_k /= 1000000
    #         # print(i, delta_k)
    #         new_k.append(sets_of_k_stiff[i]-delta_k)
    #     new_k.append(sets_of_k_stiff[8])
    #     new_k.append(sets_of_k_stiff[9])
    #     new_k.append(sets_of_k_stiff[10])
    #     if min(new_k[:8]) > 0:
    #         sets_of_k_stiff = new_k
    #     else:
    #         lr *= 0.1
    #     train_time += 1