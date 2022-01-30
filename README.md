# Введение в выч. методы

## Обновление своего репозитория:

#### 1. ссылка на удаленный репозиторий: 
 git remote add mrsu https://github.com/mrsu-ru/2022-203.git

#### 2. новая ветка (создаем и переключаемся):
 git checkout -b mrsu-main
 
#### 3. затягиваем изменения из удаленного репозитория:
 git pull mrsu main
 
#### 4. переключаемся на ветку master:
 git checkout main
 
#### 5. сливаем ветку mrsu-master в ветку master:
 git merge --no-ff mrsu-main
 
#### 6. проталкиваем изменения в удаленный репозиторий origin:
 git push origin main
 
#### 7. удаляем ветку mrsu-master:
 git branch -d mrsu-main
 
#### 8. смотрим лог изменений:
 git log --graph --decorate --all
 
 
 
