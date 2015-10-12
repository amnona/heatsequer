<?xml version="1.0" encoding="UTF-8"?>
<ui version="4.0">
 <class>Dialog</class>
 <widget class="QDialog" name="Dialog">
  <property name="geometry">
   <rect>
    <x>0</x>
    <y>0</y>
    <width>400</width>
    <height>300</height>
   </rect>
  </property>
  <property name="windowTitle">
   <string>Dialog</string>
  </property>
  <widget class="QCheckBox" name="cNumeric">
   <property name="geometry">
    <rect>
     <x>300</x>
     <y>40</y>
     <width>89</width>
     <height>20</height>
    </rect>
   </property>
   <property name="text">
    <string>Numeric</string>
   </property>
  </widget>
  <widget class="QPushButton" name="bCancel">
   <property name="geometry">
    <rect>
     <x>270</x>
     <y>220</y>
     <width>121</width>
     <height>61</height>
    </rect>
   </property>
   <property name="text">
    <string>Cancel</string>
   </property>
  </widget>
  <widget class="QLineEdit" name="tNewName">
   <property name="geometry">
    <rect>
     <x>110</x>
     <y>190</y>
     <width>161</width>
     <height>21</height>
    </rect>
   </property>
  </widget>
  <widget class="QLabel" name="label_2">
   <property name="geometry">
    <rect>
     <x>130</x>
     <y>110</y>
     <width>121</width>
     <height>71</height>
    </rect>
   </property>
   <property name="text">
    <string/>
   </property>
   <property name="pixmap">
    <pixmap>arrow_down.png</pixmap>
   </property>
   <property name="scaledContents">
    <bool>false</bool>
   </property>
  </widget>
  <widget class="QCheckBox" name="cOverwrite">
   <property name="geometry">
    <rect>
     <x>300</x>
     <y>190</y>
     <width>89</width>
     <height>20</height>
    </rect>
   </property>
   <property name="text">
    <string>OverWrite</string>
   </property>
  </widget>
  <widget class="QPushButton" name="bFieldValues">
   <property name="geometry">
    <rect>
     <x>170</x>
     <y>70</y>
     <width>61</width>
     <height>32</height>
    </rect>
   </property>
   <property name="text">
    <string>Values</string>
   </property>
  </widget>
  <widget class="QPushButton" name="bOK">
   <property name="geometry">
    <rect>
     <x>10</x>
     <y>220</y>
     <width>121</width>
     <height>61</height>
    </rect>
   </property>
   <property name="text">
    <string>Sort</string>
   </property>
  </widget>
  <widget class="QLabel" name="label">
   <property name="geometry">
    <rect>
     <x>60</x>
     <y>40</y>
     <width>41</width>
     <height>16</height>
    </rect>
   </property>
   <property name="text">
    <string>Field:</string>
   </property>
  </widget>
  <widget class="QLabel" name="lExpName">
   <property name="geometry">
    <rect>
     <x>180</x>
     <y>0</y>
     <width>59</width>
     <height>16</height>
    </rect>
   </property>
   <property name="text">
    <string>ExpName</string>
   </property>
  </widget>
  <widget class="QComboBox" name="cField">
   <property name="geometry">
    <rect>
     <x>100</x>
     <y>40</y>
     <width>191</width>
     <height>26</height>
    </rect>
   </property>
  </widget>
  <widget class="QLabel" name="lFSExpSamples_3">
   <property name="geometry">
    <rect>
     <x>40</x>
     <y>190</y>
     <width>71</width>
     <height>16</height>
    </rect>
   </property>
   <property name="text">
    <string>New name:</string>
   </property>
  </widget>
 </widget>
 <resources/>
 <connections/>
</ui>
