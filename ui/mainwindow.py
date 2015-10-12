<?xml version="1.0" encoding="UTF-8"?>
<ui version="4.0">
 <class>MainWindow</class>
 <widget class="QMainWindow" name="MainWindow">
  <property name="windowModality">
   <enum>Qt::WindowModal</enum>
  </property>
  <property name="geometry">
   <rect>
    <x>0</x>
    <y>0</y>
    <width>400</width>
    <height>300</height>
   </rect>
  </property>
  <property name="windowTitle">
   <string>MainWindow</string>
  </property>
  <widget class="QWidget" name="centralWidget">
   <widget class="QTextEdit" name="ttable">
    <property name="geometry">
     <rect>
      <x>50</x>
      <y>10</y>
      <width>201</width>
      <height>21</height>
     </rect>
    </property>
    <property name="verticalScrollBarPolicy">
     <enum>Qt::ScrollBarAlwaysOff</enum>
    </property>
    <property name="horizontalScrollBarPolicy">
     <enum>Qt::ScrollBarAlwaysOff</enum>
    </property>
    <property name="lineWrapMode">
     <enum>QTextEdit::NoWrap</enum>
    </property>
   </widget>
   <widget class="QPushButton" name="btable">
    <property name="geometry">
     <rect>
      <x>250</x>
      <y>0</y>
      <width>81</width>
      <height>41</height>
     </rect>
    </property>
    <property name="toolTip">
     <string>Choose biom table file</string>
    </property>
    <property name="text">
     <string>Choose</string>
    </property>
   </widget>
   <widget class="QLabel" name="label">
    <property name="geometry">
     <rect>
      <x>10</x>
      <y>10</y>
      <width>59</width>
      <height>16</height>
     </rect>
    </property>
    <property name="text">
     <string>Table:</string>
    </property>
   </widget>
   <widget class="QLabel" name="label_2">
    <property name="geometry">
     <rect>
      <x>10</x>
      <y>50</y>
      <width>59</width>
      <height>16</height>
     </rect>
    </property>
    <property name="text">
     <string>Table:</string>
    </property>
   </widget>
   <widget class="QPushButton" name="bmap">
    <property name="geometry">
     <rect>
      <x>250</x>
      <y>40</y>
      <width>81</width>
      <height>41</height>
     </rect>
    </property>
    <property name="toolTip">
     <string>Choose biom table file</string>
    </property>
    <property name="text">
     <string>Choose</string>
    </property>
   </widget>
   <widget class="QTextEdit" name="tmap">
    <property name="geometry">
     <rect>
      <x>50</x>
      <y>50</y>
      <width>201</width>
      <height>21</height>
     </rect>
    </property>
    <property name="verticalScrollBarPolicy">
     <enum>Qt::ScrollBarAlwaysOff</enum>
    </property>
    <property name="horizontalScrollBarPolicy">
     <enum>Qt::ScrollBarAlwaysOff</enum>
    </property>
    <property name="lineWrapMode">
     <enum>QTextEdit::NoWrap</enum>
    </property>
   </widget>
   <widget class="QPushButton" name="bload">
    <property name="geometry">
     <rect>
      <x>30</x>
      <y>180</y>
      <width>151</width>
      <height>61</height>
     </rect>
    </property>
    <property name="text">
     <string>Load</string>
    </property>
   </widget>
   <widget class="QCheckBox" name="cdeblurred">
    <property name="geometry">
     <rect>
      <x>30</x>
      <y>90</y>
      <width>89</width>
      <height>20</height>
     </rect>
    </property>
    <property name="toolTip">
     <string>biom table is from deblurring</string>
    </property>
    <property name="text">
     <string>Deblurred</string>
    </property>
    <property name="checked">
     <bool>true</bool>
    </property>
   </widget>
   <widget class="QTextEdit" name="tname">
    <property name="geometry">
     <rect>
      <x>60</x>
      <y>120</y>
      <width>201</width>
      <height>21</height>
     </rect>
    </property>
    <property name="toolTip">
     <string>Name of experiment for display (optional)</string>
    </property>
    <property name="verticalScrollBarPolicy">
     <enum>Qt::ScrollBarAlwaysOff</enum>
    </property>
    <property name="horizontalScrollBarPolicy">
     <enum>Qt::ScrollBarAlwaysOff</enum>
    </property>
    <property name="lineWrapMode">
     <enum>QTextEdit::NoWrap</enum>
    </property>
   </widget>
   <widget class="QLabel" name="label_3">
    <property name="geometry">
     <rect>
      <x>10</x>
      <y>120</y>
      <width>59</width>
      <height>16</height>
     </rect>
    </property>
    <property name="palette">
     <palette>
      <active>
       <colorrole role="WindowText">
        <brush brushstyle="SolidPattern">
         <color alpha="255">
          <red>69</red>
          <green>92</green>
          <blue>229</blue>
         </color>
        </brush>
       </colorrole>
      </active>
      <inactive>
       <colorrole role="WindowText">
        <brush brushstyle="SolidPattern">
         <color alpha="255">
          <red>69</red>
          <green>92</green>
          <blue>229</blue>
         </color>
        </brush>
       </colorrole>
      </inactive>
      <disabled>
       <colorrole role="WindowText">
        <brush brushstyle="SolidPattern">
         <color alpha="255">
          <red>127</red>
          <green>127</green>
          <blue>127</blue>
         </color>
        </brush>
       </colorrole>
      </disabled>
     </palette>
    </property>
    <property name="text">
     <string>Name:</string>
    </property>
   </widget>
   <widget class="QPushButton" name="bcancel">
    <property name="geometry">
     <rect>
      <x>220</x>
      <y>180</y>
      <width>151</width>
      <height>61</height>
     </rect>
    </property>
    <property name="text">
     <string>Cancel</string>
    </property>
   </widget>
  </widget>
  <widget class="QToolBar" name="mainToolBar">
   <attribute name="toolBarArea">
    <enum>TopToolBarArea</enum>
   </attribute>
   <attribute name="toolBarBreak">
    <bool>false</bool>
   </attribute>
  </widget>
  <widget class="QStatusBar" name="statusBar"/>
 </widget>
 <layoutdefault spacing="6" margin="11"/>
 <resources/>
 <connections/>
</ui>
